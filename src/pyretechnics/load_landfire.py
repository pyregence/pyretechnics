# [[file:../../org/pyretechnics.org::load-landfire-imports][load-landfire-imports]]
import cython
import numpy as np
import rasterio
if cython.compiled:
    from cython.cimports.pyretechnics.space_time_cube import SpaceTimeCube
else:
    from pyretechnics.space_time_cube import SpaceTimeCube
# load-landfire-imports ends here
# [[file:../../org/pyretechnics.org::load-raster][load-raster]]
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
# [[file:../../org/pyretechnics.org::load-and-convert-landfire-rasters][load-and-convert-landfire-rasters]]
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
# [[file:../../org/pyretechnics.org::verify-raster-constraints][verify-raster-constraints]]
def verify_cube_compatible_dimensions(cube_shape, rasters):
    cube_shape_ = np.asarray(cube_shape)
    for r in rasters:
        raster_shape = np.asarray((r["metadata"]["bands"],
                                   r["metadata"]["rows"],
                                   r["metadata"]["cols"]))
        if any(map(lambda x: x != 0, cube_shape_ % raster_shape)):
            raise ValueError("Some rasters do not evenly divide the space-time cube dimensions.")

    return True


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


def verify_raster_constraints(cube_shape, rasters):
    return verify_cube_compatible_dimensions(cube_shape, rasters) and verify_same_georeferences(rasters)
# verify-raster-constraints ends here
# [[file:../../org/pyretechnics.org::convert-rasters-to-space-time-cubes][convert-rasters-to-space-time-cubes]]
def convert_rasters_to_space_time_cubes(cube_shape, raster_dict):
    fn_dict = {}

    for name, raster in raster_dict.items():
        fn_dict[name] = SpaceTimeCube(cube_shape, raster["array"])

    return fn_dict
# convert-rasters-to-space-time-cubes ends here
# [[file:../../org/pyretechnics.org::read-landfire-rasters-as-space-time-cubes][read-landfire-rasters-as-space-time-cubes]]
def read_landfire_rasters_as_space_time_cubes(cube_shape, landfire_file_paths):
    landfire_rasters = load_and_convert_landfire_rasters(landfire_file_paths)
    if verify_raster_constraints(cube_shape, landfire_rasters.values()):
        return convert_rasters_to_space_time_cubes(cube_shape, landfire_rasters)
# read-landfire-rasters-as-space-time-cubes ends here
