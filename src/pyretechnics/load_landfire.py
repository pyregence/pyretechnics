# [[file:../../org/pyretechnics.org::load-landfire-imports][load-landfire-imports]]
import cython
import numpy as np
import rasterio
from rasterio.enums import Resampling
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
def maybe_resample_resolution(resolution, divisor, resampling_policy):
    if (resolution % divisor == 0):
        return resolution
    elif resampling_policy == "always_upsample":
        return divisor * (resolution // divisor + 1)
    elif resampling_policy == "nearest_match":
        return divisor * max(1, round(resolution / divisor))
    else:
        raise ValueError("The resampling_policy must be either 'always_upsample' or 'nearest_match'.")


def load_raster(file_path, dtype=None, cube_shape_divisors=(1,1,1),
                resampling_policy="nearest_match",
                resampling_method=Resampling.nearest):
    """
    resampling_policy: "always_upsample" or "nearest_match"
    resampling_method: any rasterio.enums.Resampling method
    """
    with rasterio.open(file_path, "r") as input_raster:
        metadata  = raster_metadata(input_raster)
        bands     = metadata["bands"]
        rows      = metadata["rows"]
        cols      = metadata["cols"]
        (b, r, c) = cube_shape_divisors
        new_bands = maybe_resample_resolution(bands, b, resampling_policy)
        new_rows  = maybe_resample_resolution(rows, r, resampling_policy)
        new_cols  = maybe_resample_resolution(cols, c, resampling_policy)
        if new_bands == bands and new_rows == rows and new_cols == cols:
            return {
                "array"   : input_raster.read(out_dtype=dtype),
                "metadata": metadata,
            }
        else:
            metadata["bands"]     = new_bands
            metadata["rows"]      = new_rows
            metadata["cols"]      = new_cols
            metadata["transform"] = (input_raster.transform
                                     * input_raster.transform.scale(
                                         cols / new_cols,
                                         rows / new_rows,
                                     ))
            array = input_raster.read(
                out_dtype=dtype,
                out_shape=(new_bands, new_rows, new_cols),
                resampling=resampling_method,
            )
            return {
                "array"   : array,
                "metadata": metadata,
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
        if np.any(raster_shape > cube_shape_):
            raise ValueError("Some raster dimensions exceed the space-time cube dimensions.")

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
