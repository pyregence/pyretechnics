# [[file:../../org/pyretechnics.org::test-space-time-cube][test-space-time-cube]]
import numpy as np
from pyretechnics.space_time_cube import SpaceTimeCube, LazySpaceTimeCube

#==============================================================
# Raw Data - Managed by the Caller
#==============================================================

# 2D Arrays (e.g. 30m x 30m resolution, 30km x 30km extent)
elevation_layer                    = np.arange(0,1000000).reshape(1000,1000)
slope_layer                        = np.arange(0,1000000).reshape(1000,1000)
aspect_layer                       = np.arange(0,1000000).reshape(1000,1000)
fuel_model_layer                   = np.arange(0,1000000).reshape(1000,1000)
canopy_cover_layer                 = np.arange(0,1000000).reshape(1000,1000)
canopy_height_layer                = np.arange(0,1000000).reshape(1000,1000)
canopy_base_height_layer           = np.arange(0,1000000).reshape(1000,1000)
canopy_bulk_density_layer          = np.arange(0,1000000).reshape(1000,1000)
fuel_spread_adjustment_layer       = np.arange(0,1000000).reshape(1000,1000) # Optional
suppression_difficulty_index_layer = np.arange(0,1000000).reshape(1000,1000) # Optional

# 3D Arrays (e.g. 1hr x 300m x 300m resolution, 1day x 30km x 30km extent)
temperature_layer                   = np.arange(240000).reshape(24,100,100)
relative_humidity_layer             = np.arange(240000).reshape(24,100,100)
wind_speed_10m_layer                = np.arange(240000).reshape(24,100,100)
upwind_direction_layer              = np.arange(240000).reshape(24,100,100)
fuel_moisture_dead_1hr_layer        = np.arange(240000).reshape(24,100,100)
fuel_moisture_dead_10hr_layer       = np.arange(240000).reshape(24,100,100)
fuel_moisture_dead_100hr_layer      = np.arange(240000).reshape(24,100,100)
fuel_moisture_live_herbaceous_layer = np.arange(240000).reshape(24,100,100)
fuel_moisture_live_woody_layer      = np.arange(240000).reshape(24,100,100)
foliar_moisture_layer               = np.arange(240000).reshape(24,100,100)
weather_spread_adjustment_layer     = np.arange(240000).reshape(24,100,100) # Optional

#==============================================================
# Specify the SpaceTimeCube Dimensions
#==============================================================

cube_shape = (
    24,   # bands: 1 day @ 1 hour/band
    1000, # rows:  30 km @ 30 meters/row
    1000, # cols:  30 km @ 30 meters/col
)

#==============================================================
# Create the Dictionary of Layer Names to SpaceTimeCubes
#==============================================================

def test_make_layer_lookup():
    layer_lookup = {
        # 2D Arrays (e.g. 30m x 30m resolution, 30km x 30km extent)
        "elevation"                    : SpaceTimeCube(cube_shape, elevation_layer),
        "slope"                        : SpaceTimeCube(cube_shape, slope_layer),
        "aspect"                       : SpaceTimeCube(cube_shape, aspect_layer),
        "fuel_model"                   : SpaceTimeCube(cube_shape, fuel_model_layer),
        "canopy_cover"                 : SpaceTimeCube(cube_shape, canopy_cover_layer),
        "canopy_height"                : SpaceTimeCube(cube_shape, canopy_height_layer),
        "canopy_base_height"           : SpaceTimeCube(cube_shape, canopy_base_height_layer),
        "canopy_bulk_density"          : SpaceTimeCube(cube_shape, canopy_bulk_density_layer),
        "fuel_spread_adjustment"       : SpaceTimeCube(cube_shape, fuel_spread_adjustment_layer),       # Optional
        "suppression_difficulty_index" : SpaceTimeCube(cube_shape, suppression_difficulty_index_layer), # Optional

        # 3D Arrays (e.g. 1hr x 300m x 300m resolution, 1day x 30km x 30km extent)
        "temperature"                  : SpaceTimeCube(cube_shape, temperature_layer),
        "relative_humidity"            : SpaceTimeCube(cube_shape, relative_humidity_layer),
        "wind_speed_10m"               : SpaceTimeCube(cube_shape, wind_speed_10m_layer),
        "upwind_direction"             : SpaceTimeCube(cube_shape, upwind_direction_layer),
        "fuel_moisture_dead_1hr"       : SpaceTimeCube(cube_shape, fuel_moisture_dead_1hr_layer),
        "fuel_moisture_dead_10hr"      : SpaceTimeCube(cube_shape, fuel_moisture_dead_10hr_layer),
        "fuel_moisture_dead_100hr"     : SpaceTimeCube(cube_shape, fuel_moisture_dead_100hr_layer),
        "fuel_moisture_live_herbaceous": SpaceTimeCube(cube_shape, fuel_moisture_live_herbaceous_layer),
        "fuel_moisture_live_woody"     : SpaceTimeCube(cube_shape, fuel_moisture_live_woody_layer),
        "foliar_moisture"              : SpaceTimeCube(cube_shape, foliar_moisture_layer),
        "weather_spread_adjustment"    : SpaceTimeCube(cube_shape, weather_spread_adjustment_layer),    # Optional
    }
    assert all(map(lambda cube: isinstance(cube.data, np.ndarray), layer_lookup.values()))
    return layer_lookup

#==============================================================
# Looking Up Values in the Layers
#==============================================================

def test_use_layer_lookup_2d():
    layer_lookup = test_make_layer_lookup()
    dem_100_100  = layer_lookup["elevation"].get(0,100,100)
    slp_100_100  = layer_lookup["slope"].get(0,100,100)
    asp_100_100  = layer_lookup["aspect"].get(0,100,100)
    fbfm_100_100 = layer_lookup["fuel_model"].get(0,100,100)
    cc_100_100   = layer_lookup["canopy_cover"].get(0,100,100)
    ch_100_100   = layer_lookup["canopy_height"].get(0,100,100)
    cbh_100_100  = layer_lookup["canopy_base_height"].get(0,100,100)
    cbd_100_100  = layer_lookup["canopy_bulk_density"].get(0,100,100)
    fsa_100_100  = layer_lookup["fuel_spread_adjustment"].get(0,100,100)           # Optional
    sdi_100_100  = layer_lookup["suppression_difficulty_index"].get(0,100,100)     # Optional
    assert dem_100_100  == 100100
    assert slp_100_100  == 100100
    assert asp_100_100  == 100100
    assert fbfm_100_100 == 100100
    assert cc_100_100   == 100100
    assert ch_100_100   == 100100
    assert cbh_100_100  == 100100
    assert cbd_100_100  == 100100
    assert fsa_100_100  == 100100
    assert sdi_100_100  == 100100


def test_use_layer_lookup_3d():
    layer_lookup     = test_make_layer_lookup()
    temp_12_100_100  = layer_lookup["temperature"].get(12,100,100)
    rh_12_100_100    = layer_lookup["relative_humidity"].get(12,100,100)
    wsp_12_100_100   = layer_lookup["wind_speed_10m"].get(12,100,100)
    wdir_12_100_100  = layer_lookup["upwind_direction"].get(12,100,100)
    md1_12_100_100   = layer_lookup["fuel_moisture_dead_1hr"].get(12,100,100)
    md10_12_100_100  = layer_lookup["fuel_moisture_dead_10hr"].get(12,100,100)
    md100_12_100_100 = layer_lookup["fuel_moisture_dead_100hr"].get(12,100,100)
    mlh_12_100_100   = layer_lookup["fuel_moisture_live_herbaceous"].get(12,100,100)
    mlw_12_100_100   = layer_lookup["fuel_moisture_live_woody"].get(12,100,100)
    fm_12_100_100    = layer_lookup["foliar_moisture"].get(12,100,100)
    wsa_12_100_100   = layer_lookup["weather_spread_adjustment"].get(12,100,100) # Optional
    assert temp_12_100_100  == 121010
    assert rh_12_100_100    == 121010
    assert wspx_12_100_100  == 121010
    assert wspy_12_100_100  == 121010
    assert md1_12_100_100   == 121010
    assert md10_12_100_100  == 121010
    assert md100_12_100_100 == 121010
    assert mlh_12_100_100   == 121010
    assert mlw_12_100_100   == 121010
    assert fm_12_100_100    == 121010
    assert wsa_12_100_100   == 121010
# test-space-time-cube ends here
