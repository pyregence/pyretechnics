# [[file:../../org/pyretechnics.org::add-landfire-layers-to-pyretechnics-inputs][add-landfire-layers-to-pyretechnics-inputs]]
import os
from pyretechnics.load_landfire import read_landfire_rasters_as_pyretechnics_inputs

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
    "elevation"          : project_root + "/test/data/LF2020_Elev_220_CONUS/LC20_Elev_220.tif",
    "slope"              : project_root + "/test/data/LF2020_SlpP_220_CONUS/LC20_SlpP_220.tif",
    "aspect"             : project_root + "/test/data/LF2020_Asp_220_CONUS/LC20_Asp_220.tif",
    "fuel_model"         : project_root + "/test/data/LF2022_FBFM40_230_CONUS/LC22_F40_230.tif",
    "canopy_cover"       : project_root + "/test/data/LF2022_CC_230_CONUS/LC22_CC_230.tif",
    "canopy_height"      : project_root + "/test/data/LF2022_CH_230_CONUS/LC22_CH_230.tif",
    "canopy_base_height" : project_root + "/test/data/LF2022_CBH_230_CONUS/LC22_CBH_230.tif",
    "canopy_bulk_density": project_root + "/test/data/LF2022_CBD_230_CONUS/LC22_CBD_230.tif",
}


def test_read_landfire_rasters():
    input_layer_dict = read_landfire_rasters_as_pyretechnics_inputs(landfire_file_paths)
    assert type(input_layer_dict) == dict
    assert input_layer_dict.keys() == landfire_file_paths.keys()
    assert all(map(lambda v: callable(v), input_layer_dict.values()))
    return input_layer_dict
# add-landfire-layers-to-pyretechnics-inputs ends here
# [[file:../../org/pyretechnics.org::add-constant-wind-moisture-to-pyretechnics-inputs][add-constant-wind-moisture-to-pyretechnics-inputs]]
weather_functions = {
    "wind_speed_10m_x"             : lambda t,y,x: 0.00, # meters/minute
    "wind_speed_10m_y"             : lambda t,y,x: 0.00, # meters/minute
    "fuel_moisture_dead_1hr"       : lambda t,y,x: 0.06, # ratio [0-1+] grams moisture/grams ovendry wood
    "fuel_moisture_dead_10hr"      : lambda t,y,x: 0.08, # ratio [0-1+] grams moisture/grams ovendry wood
    "fuel_moisture_dead_100hr"     : lambda t,y,x: 0.10, # ratio [0-1+] grams moisture/grams ovendry wood
    "fuel_moisture_live_herbaceous": lambda t,y,x: 0.75, # ratio [0-1+] grams moisture/grams ovendry wood
    "fuel_moisture_live_woody"     : lambda t,y,x: 0.60, # ratio [0-1+] grams moisture/grams ovendry wood
    "foliar_moisture"              : lambda t,y,x: 1.20, # ratio [0-1+] grams moisture/grams ovendry foliage
}


def test_add_weather_functions():
    input_layer_dict = test_read_landfire_rasters()
    input_layer_dict.update(weather_functions)
    assert type(input_layer_dict) == dict
    assert set(input_layer_dict.keys()) == set(landfire_file_paths.keys()).union(set(weather_functions.keys()))
    assert all(map(lambda v: callable(v), input_layer_dict.values()))
    return input_layer_dict
# add-constant-wind-moisture-to-pyretechnics-inputs ends here
# [[file:../../org/pyretechnics.org::burn-single-cell-in-pyretechnics-inputs][burn-single-cell-in-pyretechnics-inputs]]
from pyretechnics.burn_cells import compute_max_in_situ_values

def test_burn_one_cell():
    input_layer_dict = test_add_weather_functions()
    (t,y,x) = (0,100,100)
    result = compute_max_in_situ_values(input_layer_dict, t, y, x)
    assert result == {
        "max_spread_rate"        : 0.32044995422500566,
        "max_spread_direction"   : 41.0,
        "max_flame_length"       : 0.35078585296988984,
        "max_fire_line_intensity": 26.661398424207746,
        "fire_type"              : 1,
        "eccentricity"           : 0.5583790663230914,
    }
    return result
# burn-single-cell-in-pyretechnics-inputs ends here
# [[file:../../org/pyretechnics.org::burn-all-cells-in-pyretechnics-inputs][burn-all-cells-in-pyretechnics-inputs]]
import numpy as np

def test_burn_all_cells():
    # TODO: Extract these dimensions from the input layers
    rows = 613
    cols = 549

    max_spread_rate_matrix         = np.zeros((rows, cols), dtype="float32")
    max_spread_direction_matrix    = np.zeros((rows, cols), dtype="int16")
    max_flame_length_matrix        = np.zeros((rows, cols), dtype="float32")
    max_fire_line_intensity_matrix = np.zeros((rows, cols), dtype="float32")
    fire_type_matrix               = np.zeros((rows, cols), dtype="uint8")
    eccentricity_matrix            = np.zeros((rows, cols), dtype="float32")

    input_layer_dict = test_add_weather_functions()

    for y in range(rows):
        for x in range(cols):
            results                             = compute_max_in_situ_values(input_layer_dict, 0, y, x)
            max_spread_rate_matrix[y,x]         = results["max_spread_rate"]
            max_spread_direction_matrix[y,x]    = results["max_spread_direction"]
            max_flame_length_matrix[y,x]        = results["max_flame_length"]
            max_fire_line_intensity_matrix[y,x] = results["max_fire_line_intensity"]
            fire_type_matrix[y,x]               = results["fire_type"]
            eccentricity_matrix[y,x]            = results["eccentricity"]

    return {
        "max_spread_rate"        : max_spread_rate_matrix,
        "max_spread_direction"   : max_spread_direction_matrix,
        "max_flame_length"       : max_flame_length_matrix,
        "max_fire_line_intensity": max_fire_line_intensity_matrix,
        "fire_type"              : fire_type_matrix,
        "eccentricity"           : eccentricity_matrix,
    }
# burn-all-cells-in-pyretechnics-inputs ends here
