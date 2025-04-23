from argparse import ArgumentParser
import numpy as np
from pprint import pprint
import pyretechnics.eulerian_level_set as els
from pyretechnics.space_time_cube import SpaceTimeCube
import ray
import time


# NOTE: All shared arrays will have dtype=np.float32 to avoid copying
#       into the child processes when converted to SpaceTimeCubes.
def make_shared_array(shape, value):
    array     = np.full(shape, value, dtype=np.float32)
    array_ref = ray.put(array)
    return array_ref


def attach_shared_array(shared_array_ref):
    shared_array = ray.get(shared_array_ref)
    return shared_array


def empirical_LoW(fire_type_matrix):
    """
    Estimates the length/width ratio (elliptical shape) by computing the empirical covariance matrix
    of the burned coordinates and using its eigenvalues.
    """
    burned_yx = np.argwhere(fire_type_matrix > 0)
    c         = np.cov(burned_yx.astype(float), rowvar=False)
    eig       = np.linalg.eig(c)
    (v1, v2)  = eig[0] # eigenvalues
    return np.sqrt(v1 / v2)


@ray.remote
def spread_one_fire(_i, shared_array_refs, cube_shape, cube_resolution, ignited_cell,
                    start_time, max_duration, spot_config, cube_refresh_rates):
    #============================================================================================
    # Attach to the shared memory arrays specified in shared_array_refs
    #============================================================================================

    shared_arrays = {name: attach_shared_array(shared_array_ref)
                     for (name, shared_array_ref) in shared_array_refs.items()}

    #============================================================================================
    # Create an input dictionary of SpaceTimeCubes from the shared memory arrays
    #============================================================================================

    space_time_cubes = {name: SpaceTimeCube(cube_shape, shared_array)
                        for (name, shared_array) in shared_arrays.items()}

    #============================================================================================
    # Create a SpreadState object and specify a point ignition location (y, x)
    #============================================================================================

    spread_state = els.SpreadState(cube_shape).ignite_cell(ignited_cell)

    #============================================================================================
    # Spread fire from the start time for the max duration
    #============================================================================================

    runtime_start       = time.perf_counter()
    fire_spread_results = els.spread_fire_with_phi_field(space_time_cubes,
                                                         spread_state,
                                                         cube_resolution,
                                                         start_time,
                                                         max_duration,
                                                         spot_config=spot_config,
                                                         cube_refresh_rates=cube_refresh_rates)
    runtime_stop        = time.perf_counter()
    stop_time           = fire_spread_results["stop_time"]      # minutes
    stop_condition      = fire_spread_results["stop_condition"] # "max duration reached" or "no burnable cells"
    spread_state        = fire_spread_results["spread_state"]   # updated SpreadState object (mutated from inputs)
    num_tracked_cells   = fire_spread_results["num_tracked_cells"] # cell count

    #============================================================================================
    # Calculate the acres burned, total runtime, and runtime per burned cell
    #============================================================================================

    output_matrices         = spread_state.get_full_matrices(layers=["fire_type"])
    fire_type_matrix        = output_matrices["fire_type"]
    num_burned_cells        = np.count_nonzero(fire_type_matrix)             # cells
    num_crowned             = np.count_nonzero(fire_type_matrix > 1)         # cells
    acres_burned            = num_burned_cells / 4.5                         # acres
    simulation_runtime      = runtime_stop - runtime_start                   # seconds
    runtime_per_burned_cell = 1000.0 * simulation_runtime / num_burned_cells # ms/cell
    Mha_per_cell            = 0.09e-6
    ms_per_hr               = 3.6e6

    return {
        "Tracked Cells"                    : num_tracked_cells,
        "Crowning Fraction"                : num_crowned / num_burned_cells,
        "Empirical Length/Width Ratio"     : empirical_LoW(fire_type_matrix),
        "Acres Burned"                     : acres_burned,
        "Total Runtime [seconds]"          : simulation_runtime,
        "Runtime Per Burned Cell [ms/cell]": runtime_per_burned_cell,
        "Areal Throughput [Mha/(CPU.hr)]"  : f"{(ms_per_hr*Mha_per_cell/runtime_per_burned_cell):.2f}",
    }


def main(num_cores, num_jobs):
    #============================================================================================
    # Start a Ray Cluster
    #============================================================================================

    print(f"Starting Ray Cluster with {num_cores} Cores...")

    ray.init(num_cpus=num_cores)

    #============================================================================================
    # Specify the SpaceTimeCube dimensions
    #============================================================================================

    cube_shape = (
        75,   # bands: 3 days + 3 hours @ 1 hour/band
        5000, # rows:  150 km @ 30 meters/row
        5000, # cols:  150 km @ 30 meters/col
    )

    #============================================================================================
    # Specify the SpaceTimeCube resolution
    #============================================================================================

    cube_resolution = (
        60, # band_duration: minutes
        30, # cell_height:   meters
        30, # cell_width:    meters
    )

    #============================================================================================
    # Create a dictionary of shared memory arrays (using constant arrays for this example)
    #============================================================================================

    (bands, rows, cols) = cube_shape
    hi_res_grid_shape   = (rows, cols)
    lo_res_cube_shape   = (bands, rows // 10, cols // 10)

    shared_array_refs = {
        "slope"                        : make_shared_array(hi_res_grid_shape,   0.80), # rise/run
        "aspect"                       : make_shared_array(hi_res_grid_shape, 225.00), # degrees clockwise from North
        "fuel_model"                   : make_shared_array(hi_res_grid_shape, 185   ), # (TL5) high-load conifer litter
        "canopy_cover"                 : make_shared_array(hi_res_grid_shape,   0.70), # 0-1
        "canopy_height"                : make_shared_array(hi_res_grid_shape,  10.00), # m
        "canopy_base_height"           : make_shared_array(hi_res_grid_shape,   0.50), # m
        "canopy_bulk_density"          : make_shared_array(hi_res_grid_shape,   0.30), # kg/m^3
        "wind_speed_10m"               : make_shared_array(lo_res_cube_shape,  10.00), # km/hr
        "upwind_direction"             : make_shared_array(lo_res_cube_shape, 180.00), # degrees clockwise from North
        "fuel_moisture_dead_1hr"       : make_shared_array(lo_res_cube_shape,   0.03), # kg moisture/kg ovendry weight
        "fuel_moisture_dead_10hr"      : make_shared_array(lo_res_cube_shape,   0.04), # kg moisture/kg ovendry weight
        "fuel_moisture_dead_100hr"     : make_shared_array(lo_res_cube_shape,   0.05), # kg moisture/kg ovendry weight
        "fuel_moisture_live_herbaceous": make_shared_array(lo_res_cube_shape,   0.90), # kg moisture/kg ovendry weight
        "fuel_moisture_live_woody"     : make_shared_array(lo_res_cube_shape,   0.60), # kg moisture/kg ovendry weight
        "foliar_moisture"              : make_shared_array(lo_res_cube_shape,   0.70), # kg moisture/kg ovendry weight
        "temperature"                  : make_shared_array(lo_res_cube_shape,  30.00), # degrees Celsius
        "fuel_spread_adjustment"       : make_shared_array(hi_res_grid_shape,   1.00), # float >= 0.0 (Optional: default = 1.0)
        "weather_spread_adjustment"    : make_shared_array(lo_res_cube_shape,   1.00), # float >= 0.0 (Optional: default = 1.0)
    }

    #============================================================================================
    # Specify the cube refresh rates
    #============================================================================================

    cube_refresh_rates = {
        "wind_speed_10m"           : 1.0 / 15.0,
        "upwind_direction"         : 1.0 / 15.0,
        "fuel_moisture_dead_1hr"   : 1.0 / 30.0,
        "temperature"              : 1.0 / 30.0,
        "fuel_spread_adjustment"   : 0.0,
        "weather_spread_adjustment": 1.0 / 30.0,
    }

    #============================================================================================
    # Set the ignited cell, start time, and max duration of the simulation
    #============================================================================================

    # In southwestern quadrant
    ignited_cell = (500,500) # (y, x)

    # Day 2 @ 10:30am
    start_time = (24 * 60) + (10 * 60) + 30 # minutes

    # 12 hours
    max_duration = 12 * 60 # minutes

    #============================================================================================
    # Set the spotting parameters
    #============================================================================================

    spot_config = {
        "random_seed"                 : 1234567890,
        "firebrands_per_unit_heat"    : 1e-9,       # firebrands/kJ
        "downwind_distance_mean"      : 10.0,       # meters
        "fireline_intensity_exponent" : 0.3,        # downwind_distance_mean multiplier [I^fireline_intensity_exponent]
        "wind_speed_exponent"         : 0.55,       # downwind_distance_mean multiplier [U^wind_speed_exponent]
        "downwind_variance_mean_ratio": 425.0,      # meters^2 / meter [downwind_variance_mean_ratio = Var(X) / E(X)]
        "crosswind_distance_stdev"    : 100.0,      # meters
        "decay_distance"              : 200.0,      # meters
    }

    #============================================================================================
    # Collect the inputs into a tuple to be shared with the child processes
    #============================================================================================

    input_args = (
        shared_array_refs,
        cube_shape,
        cube_resolution,
        ignited_cell,
        start_time,
        max_duration,
        spot_config,
        cube_refresh_rates,
    )

    #============================================================================================
    # Run multiple fires in parallel using ray core
    #============================================================================================

    print(f"Submitting {num_jobs} Parallel Jobs...")

    futures = [spread_one_fire.remote(i, *input_args) for i in range(num_jobs)]

    results = ray.get(futures)

    for i in range(num_jobs):
        print(f"\nResults of Fire {i + 1}")
        pprint(results[i], sort_dicts=False)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--cores", type=int, help = "number of cores")
    parser.add_argument("-j", "--jobs" , type=int, help = "number of jobs")
    args = parser.parse_args()

    if args.cores and args.jobs:
        main(num_cores=args.cores, num_jobs=args.jobs)
    else:
        parser.print_help()
