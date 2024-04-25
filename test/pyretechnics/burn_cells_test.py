# [[file:../../org/Pyretechnics.org::add-landfire-layers-to-pyretechnics-inputs][add-landfire-layers-to-pyretechnics-inputs]]
from load_landfire import read_landfire_rasters_as_pyretechnics_inputs

input_layer_dict = read_landfire_rasters_as_pyretechnics_inputs({
    "elevation"          : "../data/LF2020_Elev_220_CONUS/LC20_Elev_220.tif",
    "slope"              : "../data/LF2020_SlpP_220_CONUS/LC20_SlpP_220.tif",
    "aspect"             : "../data/LF2020_Asp_220_CONUS/LC20_Asp_220.tif",
    "fuel_model"         : "../data/LF2022_FBFM40_230_CONUS/LC22_F40_230.tif",
    "canopy_cover"       : "../data/LF2022_CC_230_CONUS/LC22_CC_230.tif",
    "canopy_height"      : "../data/LF2022_CH_230_CONUS/LC22_CH_230.tif",
    "canopy_base_height" : "../data/LF2022_CBH_230_CONUS/LC22_CBH_230.tif",
    "canopy_bulk_density": "../data/LF2022_CBD_230_CONUS/LC22_CBD_230.tif",
})
# add-landfire-layers-to-pyretechnics-inputs ends here
# [[file:../../org/Pyretechnics.org::add-constant-wind-moisture-to-pyretechnics-inputs][add-constant-wind-moisture-to-pyretechnics-inputs]]
input_layer_dict.update(
    {
        "wind_speed_10m_x"             : lambda t,y,x: 0.00, # meters/minute
        "wind_speed_10m_y"             : lambda t,y,x: 0.00, # meters/minute
        "fuel_moisture_dead_1hr"       : lambda t,y,x: 0.06, # ratio [0-1+] grams moisture/grams ovendry wood
        "fuel_moisture_dead_10hr"      : lambda t,y,x: 0.08, # ratio [0-1+] grams moisture/grams ovendry wood
        "fuel_moisture_dead_100hr"     : lambda t,y,x: 0.10, # ratio [0-1+] grams moisture/grams ovendry wood
        "fuel_moisture_live_herbaceous": lambda t,y,x: 0.75, # ratio [0-1+] grams moisture/grams ovendry wood
        "fuel_moisture_live_woody"     : lambda t,y,x: 0.60, # ratio [0-1+] grams moisture/grams ovendry wood
        "foliar_moisture"              : lambda t,y,x: 1.20, # ratio [0-1+] grams moisture/grams ovendry foliage
    }
)
# add-constant-wind-moisture-to-pyretechnics-inputs ends here
# [[file:../../org/Pyretechnics.org::test-burn-cells-on-pyretechnics-inputs][test-burn-cells-on-pyretechnics-inputs]]
from burn_cells import compute_max_in_situ_values

t = 0
y = 0
x = 0

compute_max_in_situ_values(input_layer_dict, t, y, x)
# test-burn-cells-on-pyretechnics-inputs ends here
