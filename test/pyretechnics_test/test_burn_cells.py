# [[file:../../org/pyretechnics.org::add-landfire-layers-to-test-dataset][add-landfire-layers-to-test-dataset]]
import os
from pyretechnics.space_time_cube import SpaceTimeCube
from pyretechnics.load_landfire import read_landfire_rasters_as_space_time_cubes


def get_project_root(current_dir=os.curdir):
    """
    Search up the directory tree from current_dir until we find a directory containing guix.scm,
    and return this directory path. If the filesystem root directory (/) is reached, return None.
    """
    if current_dir == "/":
        return None
    elif os.path.exists(current_dir + "/guix.scm"):
        return current_dir
    else:
        return get_project_root(os.path.dirname(os.path.abspath(current_dir)))


project_root = get_project_root()


landfire_file_paths = {
    "elevation"          : project_root + "/test/data/landfire_inputs/LF2020_Elev_220_CONUS/LC20_Elev_220.tif",
    "slope"              : project_root + "/test/data/landfire_inputs/LF2020_SlpP_220_CONUS/LC20_SlpP_220.tif",
    "aspect"             : project_root + "/test/data/landfire_inputs/LF2020_Asp_220_CONUS/LC20_Asp_220.tif",
    "fuel_model"         : project_root + "/test/data/landfire_inputs/LF2022_FBFM40_230_CONUS/LC22_F40_230.tif",
    "canopy_cover"       : project_root + "/test/data/landfire_inputs/LF2022_CC_230_CONUS/LC22_CC_230.tif",
    "canopy_height"      : project_root + "/test/data/landfire_inputs/LF2022_CH_230_CONUS/LC22_CH_230.tif",
    "canopy_base_height" : project_root + "/test/data/landfire_inputs/LF2022_CBH_230_CONUS/LC22_CBH_230.tif",
    "canopy_bulk_density": project_root + "/test/data/landfire_inputs/LF2022_CBD_230_CONUS/LC22_CBD_230.tif",
}

cube_shape = (1, 613, 549) # Matches the resolution of the GeoTIFFs

def test_read_landfire_rasters():
    space_time_cubes = read_landfire_rasters_as_space_time_cubes(cube_shape, landfire_file_paths)
    assert type(space_time_cubes) == dict
    assert space_time_cubes.keys() == landfire_file_paths.keys()
    assert all(map(lambda cube: isinstance(cube, SpaceTimeCube), space_time_cubes.values()))
    return space_time_cubes
# add-landfire-layers-to-test-dataset ends here
# [[file:../../org/pyretechnics.org::add-constant-wind-and-moisture-to-test-dataset][add-constant-wind-and-moisture-to-test-dataset]]
weather_cubes = {
    "wind_speed_10m"               : SpaceTimeCube(cube_shape, 0.00), # km/hr
    "upwind_direction"             : SpaceTimeCube(cube_shape, 0.00), # degrees clockwise from North
    "fuel_moisture_dead_1hr"       : SpaceTimeCube(cube_shape, 0.06), # kg moisture/kg ovendry weight
    "fuel_moisture_dead_10hr"      : SpaceTimeCube(cube_shape, 0.08), # kg moisture/kg ovendry weight
    "fuel_moisture_dead_100hr"     : SpaceTimeCube(cube_shape, 0.10), # kg moisture/kg ovendry weight
    "fuel_moisture_live_herbaceous": SpaceTimeCube(cube_shape, 0.75), # kg moisture/kg ovendry weight
    "fuel_moisture_live_woody"     : SpaceTimeCube(cube_shape, 0.60), # kg moisture/kg ovendry weight
    "foliar_moisture"              : SpaceTimeCube(cube_shape, 1.20), # kg moisture/kg ovendry weight
}


def test_add_weather_cubes():
    space_time_cubes = test_read_landfire_rasters()
    space_time_cubes.update(weather_cubes)
    assert type(space_time_cubes) == dict
    assert set(space_time_cubes.keys()) == set(landfire_file_paths.keys()).union(set(weather_cubes.keys()))
    assert all(map(lambda cube: isinstance(cube, SpaceTimeCube), space_time_cubes.values()))
    return space_time_cubes
# add-constant-wind-and-moisture-to-test-dataset ends here
# [[file:../../org/pyretechnics.org::burn-single-cell-in-test-dataset][burn-single-cell-in-test-dataset]]
from pyretechnics.burn_cells import burn_cell_as_head_fire


def test_burn_cell_as_head_fire():
    space_time_cubes      = test_add_weather_cubes()
    space_time_coordinate = (0, 100, 100) # (t,y,x)
    spread_behavior       = burn_cell_as_head_fire(space_time_cubes,
                                                   space_time_coordinate,
                                                   surface_lw_ratio_model="rothermel")
    assert spread_behavior["fire_type"] == 1 # surface
    assert spread_behavior["spread_rate"]         - 0.32044995422500555 < 0.001
    assert spread_behavior["spread_direction"][0] - 0.644528432121562   < 0.001
    assert spread_behavior["spread_direction"][1] - 0.7414451458683358  < 0.001
    assert spread_behavior["spread_direction"][2] - 0.18666064356259804 < 0.001
    assert spread_behavior["fireline_intensity"]  - 26.66139842420774   < 0.001
    assert spread_behavior["flame_length"]        - 0.3507858529698898  < 0.001
    return spread_behavior
# burn-single-cell-in-test-dataset ends here
# [[file:../../org/pyretechnics.org::burn-all-cells-in-test-dataset][burn-all-cells-in-test-dataset]]
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.vector_utils as vu
from pyretechnics.burn_cells import burn_cell_as_head_fire


def test_burn_all_cells_as_head_fire():
    space_time_cubes     = test_add_weather_cubes()
    (_bands, rows, cols) = space_time_cubes["elevation"].shape
    grid_shape           = (rows, cols)

    max_fire_type_matrix          = np.zeros(grid_shape, dtype="uint8")
    max_spread_rate_matrix        = np.zeros(grid_shape, dtype="float32")
    max_spread_direction_matrix   = np.zeros(grid_shape, dtype="float32")
    max_fireline_intensity_matrix = np.zeros(grid_shape, dtype="float32")
    max_flame_length_matrix       = np.zeros(grid_shape, dtype="float32")

    for y in range(rows):
        for x in range(cols):
            space_time_coordinate              = (0, y, x) # (t,y,x)
            spread_behavior                    = burn_cell_as_head_fire(space_time_cubes,
                                                                        space_time_coordinate,
                                                                        surface_lw_ratio_model="rothermel")
            max_fire_type_matrix[y,x]          = spread_behavior["fire_type"]
            max_spread_rate_matrix[y,x]        = spread_behavior["spread_rate"]
            max_spread_direction_matrix[y,x]   = vu.spread_direction_vector_to_angle(spread_behavior["spread_direction"])
            max_fireline_intensity_matrix[y,x] = spread_behavior["fireline_intensity"]
            max_flame_length_matrix[y,x]       = spread_behavior["flame_length"]

    return {
        "max_fire_type"         : max_fire_type_matrix,
        "max_spread_rate"       : max_spread_rate_matrix,
        "max_spread_direction"  : max_spread_direction_matrix,
        "max_fireline_intensity": max_fireline_intensity_matrix,
        "max_flame_length"      : max_flame_length_matrix,
    }
# burn-all-cells-in-test-dataset ends here
# [[file:../../org/pyretechnics.org::load-flammap-outputs][load-flammap-outputs]]
from math import pi
from pyretechnics.load_landfire import load_raster, verify_raster_constraints


flammap_file_paths = {
    "max_fire_type"         : project_root + "/test/data/flammap_outputs/fire_type.tif",
    "max_spread_rate"       : project_root + "/test/data/flammap_outputs/ROS_ch_hr.tif",
    "max_spread_direction"  : project_root + "/test/data/flammap_outputs/max_spread_direction_radians.tif",
    "max_fireline_intensity": project_root + "/test/data/flammap_outputs/FLI_BTU_ft-s.tif",
    "max_flame_length"      : project_root + "/test/data/flammap_outputs/FL_ft.tif",
}


flammap_array_conversions = {
    #====================================================================================
    # Layer Name            : (New dtype, Mult),                # In Units -> Out Units
    #====================================================================================
    "max_fire_type"         : ("uint8"  , 1.0),                 # 0=unburned,1=surface,2=passive_crown,3=active_crown
    "max_spread_rate"       : ("float32", 0.33528),             # ch/hr -> m/min
    "max_spread_direction"  : ("float32", 180.0 / pi),          # radians -> degrees
    "max_fireline_intensity": ("float32", 3.46165186),          # Btu/ft/s -> kW/m
    "max_flame_length"      : ("float32", 0.30478512648582745), # ft -> m
}


def load_and_convert_flammap_rasters(flammap_file_paths):
    flammap_rasters = {}

    for name, path in flammap_file_paths.items():
        (dtype, multiplier)    = flammap_array_conversions[name]
        flammap_rasters[name]  = load_raster(path, dtype)
        array                  = flammap_rasters[name]["array"]
        nodata                 = flammap_rasters[name]["metadata"]["nodata"]
        array[array == nodata] = 0
        if multiplier != 1:
            array *= multiplier

    return flammap_rasters


def read_flammap_outputs(flammap_file_paths):
    cube_shape  = (1, 613, 549) # Matches the resolution of the GeoTIFFs
    raster_dict = load_and_convert_flammap_rasters(flammap_file_paths)
    if verify_raster_constraints(cube_shape, raster_dict.values()):
        return {name: raster["array"] for name, raster in raster_dict.items()}


def test_read_flammap_outputs():
    return read_flammap_outputs(flammap_file_paths)
# load-flammap-outputs ends here
