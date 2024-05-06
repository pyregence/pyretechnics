# [[file:../../org/pyretechnics.org::add-landfire-layers-to-pyretechnics-inputs][add-landfire-layers-to-pyretechnics-inputs]]
from pyretechnics.load_landfire import read_landfire_rasters_as_pyretechnics_inputs

landfire_file_paths = {
    "elevation"          : "test/data/LF2020_Elev_220_CONUS/LC20_Elev_220.tif",
    "slope"              : "test/data/LF2020_SlpP_220_CONUS/LC20_SlpP_220.tif",
    "aspect"             : "test/data/LF2020_Asp_220_CONUS/LC20_Asp_220.tif",
    "fuel_model"         : "test/data/LF2022_FBFM40_230_CONUS/LC22_F40_230.tif",
    "canopy_cover"       : "test/data/LF2022_CC_230_CONUS/LC22_CC_230.tif",
    "canopy_height"      : "test/data/LF2022_CH_230_CONUS/LC22_CH_230.tif",
    "canopy_base_height" : "test/data/LF2022_CBH_230_CONUS/LC22_CBH_230.tif",
    "canopy_bulk_density": "test/data/LF2022_CBD_230_CONUS/LC22_CBD_230.tif",
}


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


def test_read_landfire_rasters():
    input_layer_dict = read_landfire_rasters_as_pyretechnics_inputs(landfire_file_paths)
    assert type(input_layer_dict) == dict
    assert input_layer_dict.keys() == landfire_file_paths.keys()
    assert all(map(lambda v: callable(v), input_layer_dict.values()))
    input_layer_dict
# add-landfire-layers-to-pyretechnics-inputs ends here
# [[file:../../org/pyretechnics.org::add-constant-wind-moisture-to-pyretechnics-inputs][add-constant-wind-moisture-to-pyretechnics-inputs]]
def test_add_weather_functions():
    input_layer_dict = test_read_landfire_rasters()
    input_layer_dict.update(weather_functions)
    assert type(input_layer_dict) == dict
    assert set(input_layer_dict.keys()) == set(landfire_file_paths.keys()).union(set(weather_functions.keys()))
    assert all(map(lambda v: callable(v), input_layer_dict.values()))
    input_layer_dict
# add-constant-wind-moisture-to-pyretechnics-inputs ends here
# [[file:../../org/pyretechnics.org::test-burn-cells-on-pyretechnics-inputs][test-burn-cells-on-pyretechnics-inputs]]
from pyretechnics.burn_cells import compute_max_in_situ_values

def test_burn_one_cell():
    input_layer_dict = test_add_weather_functions()
    t = 0
    y = 100
    x = 100
    result = compute_max_in_situ_values(input_layer_dict, t, y, x)
    assert result == {}
# test-burn-cells-on-pyretechnics-inputs ends here
