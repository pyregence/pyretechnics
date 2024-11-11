# [[file:../org/pyretechnics.org::run-spread-fire-with-phi-field][run-spread-fire-with-phi-field]]
import numpy as np
import pyretechnics.eulerian_level_set as els
from pyretechnics.space_time_cube import SpaceTimeCube

#============================================================================================
# Specify the SpaceTimeCube dimensions
#============================================================================================

cube_shape = (
    240, # bands: 10 days @ 1 hour/band
    100, # rows:  3 km @ 30 meters/row
    100, # cols:  3 km @ 30 meters/col
)

grid_shape = cube_shape[1:]

#============================================================================================
# Specify the SpaceTimeCube resolution
#============================================================================================

cube_resolution = (
    60, # band_duration: minutes
    30, # cell_height:   meters
    30, # cell_width:    meters
)

#============================================================================================
# Create an input dictionary of SpaceTimeCubes (using constant data for this example)
#============================================================================================

space_time_cubes = {
    "slope"                        : SpaceTimeCube(cube_shape, 0.8),   # rise/run
    "aspect"                       : SpaceTimeCube(cube_shape, 225.0), # degrees clockwise from North
    "fuel_model"                   : SpaceTimeCube(cube_shape, 101),   # integer index in fm.fuel_model_table
    "canopy_cover"                 : SpaceTimeCube(cube_shape, 0.6),   # 0-1
    "canopy_height"                : SpaceTimeCube(cube_shape, 30.0),  # m
    "canopy_base_height"           : SpaceTimeCube(cube_shape, 3.0),   # m
    "canopy_bulk_density"          : SpaceTimeCube(cube_shape, 0.3),   # kg/m^3
    "wind_speed_10m"               : SpaceTimeCube(cube_shape, 30.0),  # km/hr
    "upwind_direction"             : SpaceTimeCube(cube_shape, 180.0), # degrees clockwise from North
    "fuel_moisture_dead_1hr"       : SpaceTimeCube(cube_shape, 0.05),  # kg moisture/kg ovendry weight
    "fuel_moisture_dead_10hr"      : SpaceTimeCube(cube_shape, 0.10),  # kg moisture/kg ovendry weight
    "fuel_moisture_dead_100hr"     : SpaceTimeCube(cube_shape, 0.15),  # kg moisture/kg ovendry weight
    "fuel_moisture_live_herbaceous": SpaceTimeCube(cube_shape, 0.90),  # kg moisture/kg ovendry weight
    "fuel_moisture_live_woody"     : SpaceTimeCube(cube_shape, 0.60),  # kg moisture/kg ovendry weight
    "foliar_moisture"              : SpaceTimeCube(cube_shape, 0.90),  # kg moisture/kg ovendry weight
    "fuel_spread_adjustment"       : SpaceTimeCube(cube_shape, 1.0),   # float >= 0.0 (Optional: defaults to 1.0)
    "weather_spread_adjustment"    : SpaceTimeCube(cube_shape, 1.0),   # float >= 0.0 (Optional: defaults to 1.0)
}

#============================================================================================
# Create an output dictionary of 2D Numpy arrays
#============================================================================================

output_matrices = {
    "phi"               : np.ones(grid_shape, dtype="float32"),       # 2D float array of values in [-1,1]
    "fire_type"         : np.zeros(grid_shape, dtype="uint8"),        # 2D byte array (0-3)
    "spread_rate"       : np.zeros(grid_shape, dtype="float32"),      # 2D float array (m/min)
    "spread_direction"  : np.zeros(grid_shape, dtype="float32"),      # 2D float array (degrees clockwise from North)
    "fireline_intensity": np.zeros(grid_shape, dtype="float32"),      # 2D float array (kW/m)
    "flame_length"      : np.zeros(grid_shape, dtype="float32"),      # 2D float array (m)
    "time_of_arrival"   : np.full(grid_shape, -1.0, dtype="float32"), # 2D float array (min)
}

#============================================================================================
# Set the start time, max duration, and initially ignited cell
#============================================================================================

# Day 2 @ 10:30am
start_time = 2070  # minutes

# 8 hours
max_duration = 480 # minutes

# Burn initially ignited cell into the phi matrix by setting it to -1.0
output_matrices["phi"][50,50] = -1.0

#============================================================================================
# Spread fire from the start time for the max duration
#============================================================================================

fire_spread_results = els.spread_fire_with_phi_field(space_time_cubes, output_matrices, cube_resolution,
                                                     start_time, max_duration)

stop_time       = fire_spread_results["stop_time"]       # minutes
stop_condition  = fire_spread_results["stop_condition"]  # "max duration reached" or "no burnable cells"
output_matrices = fire_spread_results["output_matrices"] # updated 2D arrays (mutated from inputs)
# run-spread-fire-with-phi-field ends here
