# [[file:../../org/pyretechnics.org::lazy-array-lookup-functions][lazy-array-lookup-functions]]
import numpy as np

#==============================================================
# Lazy-Loading Multi-Resolution Array Lookup Functions
#==============================================================

def translate_coords_2d(coords_2d, simulation_shape_2d, layer_shape_2d, chunk_shape_2d):
    """
    Given row y and column x on the simulation grid, return chunk
    row cy, chunk column cx, subchunk row c_y, and subchunk column c_x on
    the chunk grid.
    NOTE: (y,x) = (0,0) is the upper-left corner of the array.
    NOTE: simulation_shape_2d >= layer_shape_2d >= chunk_shape_2d
    """
    (y, x) = coords_2d
    (simulation_rows, simulation_cols) = simulation_shape_2d
    (layer_rows, layer_cols) = layer_shape_2d
    (chunk_rows, chunk_cols) = chunk_shape_2d
    ly = y // (simulation_rows // layer_rows)
    lx = x // (simulation_cols // layer_cols)
    cy  = ly // chunk_rows
    cx  = lx // chunk_cols
    c_y = ly % chunk_rows
    c_x = lx % chunk_cols
    return (cy, cx, c_y, c_x)

def translate_coords_3d(coords_3d, simulation_shape_3d, layer_shape_3d, chunk_shape_3d):
    """
    Given band z, row y, and column x on the simulation grid,
    return chunk band cz, chunk row cy, chunk column cx, subchunk band
    c_z, subchunk row c_y, and subchunk column c_x on the chunk grid.
    NOTE: (z,y,x) = (0,0,0) is the upper-left corner of the array in the first timestep.
    NOTE: simulation_shape_3d >= layer_shape_3d >= chunk_shape_3d
    """
    (z, y, x) = coords_3d
    (simulation_bands, simulation_rows, simulation_cols) = simulation_shape_3d
    (layer_bands, layer_rows, layer_cols) = layer_shape_3d
    (chunk_bands, chunk_rows, chunk_cols) = chunk_shape_3d
    lz = z // (simulation_bands // layer_bands)
    ly = y // (simulation_rows // layer_rows)
    lx = x // (simulation_cols // layer_cols)
    cz  = lz // chunk_bands
    cy  = ly // chunk_rows
    cx  = lx // chunk_cols
    c_z = lz % chunk_bands
    c_y = ly % chunk_rows
    c_x = lx % chunk_cols
    return (cz, cy, cx, c_z, c_y, c_x)

def get_or_load_chunk_2d(chunk_cache_2d, cy, cx, load_chunk_2d):
    """
    Checks whether chunk_cache_2d contains chunk (cy, cx). If so,
    return it. Otherwise, call load_chunk_2d(cy, cx) to retrieve it,
    store it in chunk_cache_2d, and return it.
    """
    chunk_2d = chunk_cache_2d[cy, cx]
    if (type(chunk_2d) == np.ndarray):
        return chunk_2d
    else:
        chunk_2d = load_chunk_2d(cy, cx)
        chunk_cache_2d[cy, cx] = chunk_2d
        return chunk_2d

def get_or_load_chunk_3d(chunk_cache_3d, cz, cy, cx, load_chunk_3d):
    """
    Checks whether chunk_cache_3d contains chunk (cz, cy, cx). If so,
    return it. Otherwise, call load_chunk_3d(cz, cy, cx) to retrieve it,
    store it in chunk_cache_3d, and return it.
    """
    chunk_3d = chunk_cache_3d[cz, cy, cx]
    if (type(chunk_3d) == np.ndarray):
        return chunk_3d
    else:
        chunk_3d = load_chunk_3d(cz, cy, cx)
        chunk_cache_3d[cz, cy, cx] = chunk_3d
        return chunk_3d

def at_coords_2d(coords_2d, simulation_shape_2d, layer_shape_2d, chunk_shape_2d, chunk_cache_2d, load_chunk_2d):
    """
    Given row y and column x on the simulation grid, return the
    value at that index in the underlying chunk cache.
    NOTE: (y,x) = (0,0) is the upper-left corner of the array.
    """
    (cy, cx, c_y, c_x) = translate_coords_2d(coords_2d, simulation_shape_2d, layer_shape_2d, chunk_shape_2d)
    chunk_2d = get_or_load_chunk_2d(chunk_cache_2d, cy, cx, load_chunk_2d)
    return chunk_2d[c_y, c_x]

def at_coords_3d(coords_3d, simulation_shape_3d, layer_shape_3d, chunk_shape_3d, chunk_cache_3d, load_chunk_3d):
    """
    Given band z, row y, and column x on the simulation grid,
    return the value at that index in the underlying chunk cache.
    NOTE: (z,y,x) = (0,0,0) is the upper-left corner of the array in the first timestep.
    """
    (cz, cy, cx, c_z, c_y, c_x) = translate_coords_3d(coords_3d, simulation_shape_3d, layer_shape_3d, chunk_shape_3d)
    chunk_3d = get_or_load_chunk_3d(chunk_cache_3d, cz, cy, cx, load_chunk_3d)
    return chunk_3d[c_z, c_y, c_x]

#==============================================================
# Constructor: Returns a 2D/3D Lookup Function for One Array
#==============================================================

def make_lookup_fn_2d(simulation_shape_2d, layer_shape_2d, chunk_shape_2d, load_chunk_2d):
    """
    Given the array shapes of the simulation space (|Y|,|X|), the
    underlying data layer (|LY|,|LX|), and a single chunk within the
    chunk cache (|CY|,|CX|) as well as a function to load one chunk
    on demand, return a closure that will retrieve the value from the
    underlying data layer corresponding to coordinate (y,x) in the
    simulation space. Chunks will be loaded on demand using load_chunk_2d.
    NOTE: (y,x) = (0,0) is the upper-left corner of the array.
    NOTE: simulation_shape_2d >= layer_shape_2d >= chunk_shape_2d
    """
    (layer_rows, layer_cols) = layer_shape_2d
    (chunk_rows, chunk_cols) = chunk_shape_2d
    chunk_cache_2d = np.empty((layer_rows // chunk_rows,
                               layer_cols // chunk_cols),
                              dtype=object)
    return lambda y, x: at_coords_2d((y, x),
                                     simulation_shape_2d,
                                     layer_shape_2d,
                                     chunk_shape_2d,
                                     chunk_cache_2d,
                                     load_chunk_2d)

def make_lookup_fn_3d(simulation_shape_3d, layer_shape_3d, chunk_shape_3d, load_chunk_3d):
    """
    Given the array shapes of the simulation space (|Z|,|Y|,|X|), the
    underlying data layer (|LZ|,|LY|,|LX|), and a single chunk within the
    chunk cache (|CZ|,|CY|,|CX|) as well as a function to load one chunk
    on demand, return a closure that will retrieve the value from the
    underlying data layer corresponding to coordinate (z,y,x) in the
    simulation space. Chunks will be loaded on demand using load_chunk_3d.
    NOTE: (z,y,x) = (0,0,0) is the upper-left corner of the array in the first timestep.
    NOTE: simulation_shape_3d >= layer_shape_3d >= chunk_shape_3d
    """
    (layer_bands, layer_rows, layer_cols) = layer_shape_3d
    (chunk_bands, chunk_rows, chunk_cols) = chunk_shape_3d
    chunk_cache_3d = np.empty((layer_bands // chunk_bands,
                               layer_rows // chunk_rows,
                               layer_cols // chunk_cols),
                              dtype=object)
    return lambda z, y, x: at_coords_3d((z, y, x),
                                        simulation_shape_3d,
                                        layer_shape_3d,
                                        chunk_shape_3d,
                                        chunk_cache_3d,
                                        load_chunk_3d)
# lazy-array-lookup-functions ends here
