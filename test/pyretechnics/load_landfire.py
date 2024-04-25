# [[file:../../org/Pyretechnics.org::load-raster][load-raster]]
import rasterio

def raster_metadata(raster):
    return {
        "name"      : raster.name,
        "bands"     : raster.count,
        "rows"      : raster.height,
        "cols"      : raster.width,
        "band_types": {i: dtype for i, dtype in zip(raster.indexes, raster.dtypes)},
        "nodata"    : raster.nodata,
        "crs"       : raster.crs,
        "bounds"    : raster.bounds,
        "transform" : raster.transform,
    }


# TODO: rasterio.open can also load chunks of a raster file: https://pypi.org/project/rasterio/
def load_raster(file_path, dtype=None):
    with rasterio.open(file_path, "r") as input_raster:
        return {
            "array"   : input_raster.read(out_dtype=dtype),
            "metadata": raster_metadata(input_raster),
        }
# load-raster ends here
# [[file:../../org/Pyretechnics.org::load-and-convert-landfire-rasters][load-and-convert-landfire-rasters]]
landfire_array_conversions = {
    #==============================================================================
    # Layer Name         : (New dtype, Mult), # New Units                [Min-Max]
    #==============================================================================
    "elevation"          : ("float32", 1.00), # meters above sea level   [0-8850]
    "slope"              : ("float32", 0.01), # meters rise / meters run [0-4]
    "aspect"             : ("int16"  , 1   ), # degrees CW from North    [0-359]
    "fuel_model"         : ("int16"  , 1   ), # integer classes          [1-204]
    "canopy_cover"       : ("float32", 0.01), # canopy area / total area [0-0.95]
    "canopy_height"      : ("float32", 0.10), # meters                   [0-51]
    "canopy_base_height" : ("float32", 0.10), # meters                   [0-10]
    "canopy_bulk_density": ("float32", 0.01), # kilograms/meters^3       [0-0.45]
}


def load_and_convert_landfire_rasters(landfire_file_paths):
    landfire_rasters = {}

    for name, path in landfire_file_paths.items():
        (dtype, multiplier) = landfire_array_conversions[name]
        landfire_rasters[name] = load_raster(path, dtype)
        if multiplier != 1:
            array  = landfire_rasters[name]["array"]
            nodata = landfire_rasters[name]["metadata"]["nodata"]
            array[array != nodata] *= multiplier

    return landfire_rasters
# load-and-convert-landfire-rasters ends here
# [[file:../../org/Pyretechnics.org::verify-raster-constraints][verify-raster-constraints]]
def verify_same_dimensions(rasters):
    dimensions = [
        (r["metadata"]["bands"],
         r["metadata"]["rows"],
         r["metadata"]["cols"])
        for r in rasters
    ]
    if len(set(dimensions)) == 1:
        return True
    else:
        raise ValueError("All rasters do not share the same dimensions.")


def verify_same_georeferences(rasters):
    georeferences = [
        (r["metadata"]["crs"],
         r["metadata"]["bounds"],
         r["metadata"]["transform"])
        for r in rasters
    ]
    if len(set(georeferences)) == 1:
        return True
    else:
        raise ValueError("All rasters do not share the same georeferences.")


def verify_raster_constraints(rasters):
    return verify_same_dimensions(rasters) and verify_same_georeferences(rasters)
# verify-raster-constraints ends here
# [[file:../../org/Pyretechnics.org::convert-rasters-to-2d-fns][convert-rasters-to-2d-fns]]
# TODO: This does not use the lazy chunking, multi-resolution functions defined earlier.
def convert_rasters_to_2d_fns(raster_dict):
    fn_dict = {}

    for name, raster in raster_dict.items():
        fn_dict[name] = (lambda array: lambda t,y,x: array[0,y,x])(raster["array"])

    return fn_dict
# convert-rasters-to-2d-fns ends here
# [[file:../../org/Pyretechnics.org::read-landfire-rasters-as-pyretechnics-inputs][read-landfire-rasters-as-pyretechnics-inputs]]
def read_landfire_rasters_as_pyretechnics_inputs(landfire_file_paths):
    landfire_rasters = load_and_convert_landfire_rasters(landfire_file_paths)
    if verify_raster_constraints(landfire_rasters.values()):
        return convert_rasters_to_2d_fns(landfire_rasters)


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
# read-landfire-rasters-as-pyretechnics-inputs ends here
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
