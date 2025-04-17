import numpy as np
import pyretechnics.eulerian_level_set as els
from pyretechnics.space_time_cube import SpaceTimeCube
import time

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
# Create an input dictionary of SpaceTimeCubes (using constant arrays for this example)
#============================================================================================

def arr2d(value):
    grid_shape = cube_shape[1:]
    return np.full(grid_shape, value, dtype=np.float32)

def arr3d(value):
    (b, r, c) = cube_shape
    return np.full((b, r//10, c//10), value, dtype=np.float32)

space_time_cubes = {
    "slope"                        : SpaceTimeCube(cube_shape, arr2d(0.8)),   # rise/run
    "aspect"                       : SpaceTimeCube(cube_shape, arr2d(225.0)), # degrees clockwise from North
    # Fuel model 185: (TL5) high-load conifer litter
    "fuel_model"                   : SpaceTimeCube(cube_shape, arr2d(185)),   # integer index in fm.fuel_model_table
    "canopy_cover"                 : SpaceTimeCube(cube_shape, arr2d(0.7)),   # 0-1
    "canopy_height"                : SpaceTimeCube(cube_shape, arr2d(10.0)),  # m
    "canopy_base_height"           : SpaceTimeCube(cube_shape, arr2d(0.5)),   # m
    "canopy_bulk_density"          : SpaceTimeCube(cube_shape, arr2d(0.3)),   # kg/m^3
    "wind_speed_10m"               : SpaceTimeCube(cube_shape, arr3d(10.0)),  # km/hr
    "upwind_direction"             : SpaceTimeCube(cube_shape, arr3d(180.0)), # degrees clockwise from North
    "fuel_moisture_dead_1hr"       : SpaceTimeCube(cube_shape, arr3d(0.03)),  # kg moisture/kg ovendry weight
    "fuel_moisture_dead_10hr"      : SpaceTimeCube(cube_shape, arr3d(0.04)),  # kg moisture/kg ovendry weight
    "fuel_moisture_dead_100hr"     : SpaceTimeCube(cube_shape, arr3d(0.05)),  # kg moisture/kg ovendry weight
    "fuel_moisture_live_herbaceous": SpaceTimeCube(cube_shape, arr3d(0.90)),  # kg moisture/kg ovendry weight
    "fuel_moisture_live_woody"     : SpaceTimeCube(cube_shape, arr3d(0.60)),  # kg moisture/kg ovendry weight
    "foliar_moisture"              : SpaceTimeCube(cube_shape, arr3d(0.70)),  # kg moisture/kg ovendry weight
    "temperature"                  : SpaceTimeCube(cube_shape, arr3d(30.0)),  # degrees Celsius
    "fuel_spread_adjustment"       : SpaceTimeCube(cube_shape, arr2d(1.0)),   # float >= 0.0 (Optional: default = 1.0)
    "weather_spread_adjustment"    : SpaceTimeCube(cube_shape, arr3d(1.0)),   # float >= 0.0 (Optional: default = 1.0)
}

cube_refresh_rates = {
    "wind_speed_10m"           : 1.0 / 15.0,
    "upwind_direction"         : 1.0 / 15.0,
    "fuel_moisture_dead_1hr"   : 1.0 / 30.0,
    "temperature"              : 1.0 / 30.0,
    "fuel_spread_adjustment"   : 0.0,
    "weather_spread_adjustment": 1.0 / 30.0,
}

#============================================================================================
# Create a SpreadState object and specify a point ignition location (y, x)
#============================================================================================

spread_state = els.SpreadState(cube_shape).ignite_cell((500,500))

#============================================================================================
# Set the start time and max duration of the simulation
#============================================================================================

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
stop_time           = fire_spread_results["stop_time"]         # minutes
stop_condition      = fire_spread_results["stop_condition"]    # "max duration reached" or "no burnable cells"
spread_state        = fire_spread_results["spread_state"]      # updated SpreadState object (mutated from inputs)
num_tracked_cells   = fire_spread_results["num_tracked_cells"] # cell count

#============================================================================================
# Print out the acres burned, total runtime, and runtime per burned cell
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

def empirical_LoW(fire_type_matrix):
    """
    Estimates the length/width ratio (elliptical shape) by computing the empirical covariance matrix
    of the burned coordinates and using its eigenvalues.
    """
    burned_yx = np.argwhere(fire_type_matrix > 0)
    c         = np.cov(burned_yx.astype(float), rowvar=False)
    eig       = np.linalg.eig(c)
    v1, v2    = eig[0] # eigenvalues
    return np.sqrt(v1 / v2)

print(f"Tracked Cells: {num_tracked_cells}")
print(f"Crowning Fraction: {num_crowned / num_burned_cells}")
print(f"Empirical Length/Width Ratio: {empirical_LoW(fire_type_matrix)}")
print("Acres Burned: " + str(acres_burned))
print("Total Runtime: " + str(simulation_runtime) + " seconds")
print("Runtime Per Burned Cell: " + str(runtime_per_burned_cell) + " ms/cell")
print(f"Areal Throughput: {(ms_per_hr*Mha_per_cell/runtime_per_burned_cell):.2f} Mha/(CPU.hr)")
