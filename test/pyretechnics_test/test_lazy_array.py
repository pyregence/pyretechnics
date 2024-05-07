# [[file:../../org/pyretechnics.org::*Lazy Array Usage Examples][Lazy Array Usage Examples:1]]
import numpy as np
from pyretechnics.lazy_array import make_lookup_fn_2d, make_lookup_fn_3d

#==============================================================
# Example Chunk Loading Functions
#==============================================================

def make_load_chunk_2d(layer_2d, chunk_shape_2d):
    """
    Example 2D chunk loading function. Since we are calling the
    lazy array functions from Python and they will be used as
    read-only data by the caller, we can just provide a sliced view of
    the underlying array rather than generating a copy.
    """
    (chunk_rows, chunk_cols) = chunk_shape_2d
    return lambda cy, cx: layer_2d[cy * chunk_rows : (cy + 1) * chunk_rows,
                                   cx * chunk_cols : (cx + 1) * chunk_cols]


def make_load_chunk_3d(layer_3d, chunk_shape_3d):
    """
    Example 3D chunk loading function. Since we are calling the
    lazy array functions from Python and they will be used as
    read-only data by the caller, we can just provide a sliced view of
    the underlying array rather than generating a copy.
    """
    (chunk_bands, chunk_rows, chunk_cols) = chunk_shape_3d
    return lambda cz, cy, cx: layer_3d[cz * chunk_bands : (cz + 1) * chunk_bands,
                                       cy * chunk_rows  : (cy + 1) * chunk_rows,
                                       cx * chunk_cols  : (cx + 1) * chunk_cols]

#==============================================================
# Example Lookup Functions for 1000x1000 and 24x100x100 Layers
#==============================================================

#--------------------  bands, rows, cols
simulation_shape_2d = (       1000, 1000)
layer_shape_2d      = (       1000, 1000)
chunk_shape_2d      = (        100,  100)

simulation_shape_3d = (   24, 1000, 1000)
layer_shape_3d      = (   24,  100,  100)
chunk_shape_3d      = (    1,   10,   10)

# Partial Application Functions
def make_lookup_fn_2d_for_layer(layer_2d):
    return make_lookup_fn_2d(simulation_shape_2d,
                             layer_shape_2d,
                             chunk_shape_2d,
                             make_load_chunk_2d(layer_2d, chunk_shape_2d))


def make_lookup_fn_3d_for_layer(layer_3d):
    return make_lookup_fn_3d(simulation_shape_3d,
                             layer_shape_3d,
                             chunk_shape_3d,
                             make_load_chunk_3d(layer_3d, chunk_shape_3d))

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
wind_speed_10m_x_layer              = np.arange(240000).reshape(24,100,100)
wind_speed_10m_y_layer              = np.arange(240000).reshape(24,100,100)
fuel_moisture_dead_1hr_layer        = np.arange(240000).reshape(24,100,100)
fuel_moisture_dead_10hr_layer       = np.arange(240000).reshape(24,100,100)
fuel_moisture_dead_100hr_layer      = np.arange(240000).reshape(24,100,100)
fuel_moisture_live_herbaceous_layer = np.arange(240000).reshape(24,100,100)
fuel_moisture_live_woody_layer      = np.arange(240000).reshape(24,100,100)
foliar_moisture_layer               = np.arange(240000).reshape(24,100,100)
weather_spread_adjustment_layer     = np.arange(240000).reshape(24,100,100) # Optional

#==============================================================
# Creating the Dictionary of Layer Names to Lookup Functions
#==============================================================

def test_make_layer_lookup():
    layer_lookup = {
        # 2D Arrays (e.g. 30m x 30m resolution, 30km x 30km extent)
        "elevation"                    : make_lookup_fn_2d_for_layer(elevation_layer),
        "slope"                        : make_lookup_fn_2d_for_layer(slope_layer),
        "aspect"                       : make_lookup_fn_2d_for_layer(aspect_layer),
        "fuel_model"                   : make_lookup_fn_2d_for_layer(fuel_model_layer),
        "canopy_cover"                 : make_lookup_fn_2d_for_layer(canopy_cover_layer),
        "canopy_height"                : make_lookup_fn_2d_for_layer(canopy_height_layer),
        "canopy_base_height"           : make_lookup_fn_2d_for_layer(canopy_base_height_layer),
        "canopy_bulk_density"          : make_lookup_fn_2d_for_layer(canopy_bulk_density_layer),
        "fuel_spread_adjustment"       : make_lookup_fn_2d_for_layer(fuel_spread_adjustment_layer),       # Optional
        "suppression_difficulty_index" : make_lookup_fn_2d_for_layer(suppression_difficulty_index_layer), # Optional

        # 3D Arrays (e.g. 1hr x 300m x 300m resolution, 1day x 30km x 30km extent)
        "temperature"                  : make_lookup_fn_3d_for_layer(temperature_layer),
        "relative_humidity"            : make_lookup_fn_3d_for_layer(relative_humidity_layer),
        "wind_speed_10m_x"             : make_lookup_fn_3d_for_layer(wind_speed_10m_x_layer),
        "wind_speed_10m_y"             : make_lookup_fn_3d_for_layer(wind_speed_10m_y_layer),
        "fuel_moisture_dead_1hr"       : make_lookup_fn_3d_for_layer(fuel_moisture_dead_1hr_layer),
        "fuel_moisture_dead_10hr"      : make_lookup_fn_3d_for_layer(fuel_moisture_dead_10hr_layer),
        "fuel_moisture_dead_100hr"     : make_lookup_fn_3d_for_layer(fuel_moisture_dead_100hr_layer),
        "fuel_moisture_live_herbaceous": make_lookup_fn_3d_for_layer(fuel_moisture_live_herbaceous_layer),
        "fuel_moisture_live_woody"     : make_lookup_fn_3d_for_layer(fuel_moisture_live_woody_layer),
        "foliar_moisture"              : make_lookup_fn_3d_for_layer(foliar_moisture_layer),
        "weather_spread_adjustment"    : make_lookup_fn_3d_for_layer(weather_spread_adjustment_layer),    # Optional
    }
    assert all(map(lambda v: callable(v), layer_lookup.values()))
    return layer_lookup

#==============================================================
# Looking Up Values in the Layers
#==============================================================

# NOTE: 2D coords should be provided as (y,x) in simulation space.
def test_use_layer_lookup_2d():
    layer_lookup = test_make_layer_lookup()
    dem_100_100  = layer_lookup["elevation"](100,100)
    slp_100_100  = layer_lookup["slope"](100,100)
    asp_100_100  = layer_lookup["aspect"](100,100)
    fbfm_100_100 = layer_lookup["fuel_model"](100,100)
    cc_100_100   = layer_lookup["canopy_cover"](100,100)
    ch_100_100   = layer_lookup["canopy_height"](100,100)
    cbh_100_100  = layer_lookup["canopy_base_height"](100,100)
    cbd_100_100  = layer_lookup["canopy_bulk_density"](100,100)
    fsa_100_100  = layer_lookup["fuel_spread_adjustment"](100,100)           # Optional
    sdi_100_100  = layer_lookup["suppression_difficulty_index"](100,100)     # Optional
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


# NOTE: 3D coords should be provided as (t,y,x) in simulation space.
def test_use_layer_lookup_3d():
    layer_lookup = test_make_layer_lookup()
    temp_12_100_100  = layer_lookup["temperature"](12,100,100)
    rh_12_100_100    = layer_lookup["relative_humidity"](12,100,100)
    wspx_12_100_100  = layer_lookup["wind_speed_10m_x"](12,100,100)
    wspy_12_100_100  = layer_lookup["wind_speed_10m_y"](12,100,100)
    md1_12_100_100   = layer_lookup["fuel_moisture_dead_1hr"](12,100,100)
    md10_12_100_100  = layer_lookup["fuel_moisture_dead_10hr"](12,100,100)
    md100_12_100_100 = layer_lookup["fuel_moisture_dead_100hr"](12,100,100)
    mlh_12_100_100   = layer_lookup["fuel_moisture_live_herbaceous"](12,100,100)
    mlw_12_100_100   = layer_lookup["fuel_moisture_live_woody"](12,100,100)
    fm_12_100_100    = layer_lookup["foliar_moisture"](12,100,100)
    wsa_12_100_100   = layer_lookup["weather_spread_adjustment"](12,100,100) # Optional
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
# Lazy Array Usage Examples:1 ends here
