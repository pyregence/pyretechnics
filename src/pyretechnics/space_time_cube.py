# [[file:../../org/pyretechnics.org::space-time-cube-imports][space-time-cube-imports]]
from functools import reduce
import numpy as np
# space-time-cube-imports ends here
# [[file:../../org/pyretechnics.org::space-time-cube-utilities][space-time-cube-utilities]]
def is_pos_int(x):
    return isinstance(x, int) and x > 0


def divide_evenly(dividend, divisor):
    (quotient, remainder) = divmod(dividend, divisor)
    if remainder == 0:
        return quotient
    else:
        raise ValueError(str(dividend) + " must be an exact multiple of " + str(divisor) + ".")


def to_positive_index_range(index_range, axis_length):
    """
    Translate None and negative indices to positive indices.
    """
    if index_range == None:
        return (0, axis_length)
    else:
        (start, stop) = index_range
        return (
            0 if start == None else axis_length + start if start < 0 else start,
            axis_length if stop == None else axis_length + stop if stop < 0 else stop
        )


def maybe_repeat_array(array, axis_repetitions):
    """
    Return a new array that is created by repeating the elements from the input
    array repetitions times along the specified array axis. Avoid allocating
    new memory if repetitions == 1 or if the repeated array axis has length 1.
    """
    (axis, repetitions) = axis_repetitions
    if repetitions == 1:
        return array
    else:
        array_shape = list(np.shape(array))
        if array_shape[axis] == 1:
            array_shape[axis] = repetitions
            return np.broadcast_to(array, array_shape)
        else:
            return np.repeat(array, repetitions, axis)
# space-time-cube-utilities ends here
# [[file:../../org/pyretechnics.org::space-time-cube-class][space-time-cube-class]]
import cython
if cython.compiled:
    from cython.cimports.pyretechnics.cy_types import pyidx
else:
    from pyretechnics.py_types import pyidx
import cython as cy

class ISpaceTimeCube:
    def get(self, t, x, y):
        pass
    def get_float(self, t, x, y):
        pass

@cy.cclass
class SpaceTimeCube(ISpaceTimeCube):
    """
    Create an object that represents a 3D array with dimensions (T,Y,X) given by cube_shape.
    Internally, data is stored as a 3D Numpy array at the resolution of the provided base data.
    Whenever a point value or contiguous space-time region of values is requested, translate
    the given cube_shape coordinates into base coordinates, look up the values from the base data,
    expand them (if necessary) back into the cube_shape resolution, and return the resulting scalar
    value or array to the caller.
    """
    ndim: cy.int
    size: cy.ulong
    shape: object
    base: object
    t_repetitions: cy.int
    y_repetitions: cy.int
    x_repetitions: cy.int
    data: cy.float[:,:,:]

    def __init__(self, cube_shape, base):
        """
        NOTE: The resolutions in cube_shape must be exact multiples of any existing dimensions
              in the base data.
        """
        # Ensure that cube_shape contains 3 values or throw an error
        (cube_bands, cube_rows, cube_cols) = cube_shape

        # Ensure that cube_shape only contains positive integers or throw an error
        if not(all(map(is_pos_int, cube_shape))):
            raise ValueError("The cube_shape must only contain positive integers.")


        base = np.asarray(base, dtype=np.float32)

        # Store the cube metadata for later
        self.ndim  = 3
        self.size  = cube_bands * cube_rows * cube_cols
        self.shape = cube_shape
        self.base  = base

        # Store the base data as a 3D array along with its axis repetitions
        base_dimensions = np.ndim(base)

        data: cy.float[:,:,:]
        if base_dimensions == 0:
            # 0D: Constant Input
            self.t_repetitions = cube_bands
            self.y_repetitions = cube_rows
            self.x_repetitions = cube_cols
            data =  np.asarray([[[base]]])

        elif base_dimensions == 1:
            # 1D: Time-Series Input
            base_bands = len(base)
            self.t_repetitions = divide_evenly(cube_bands, base_bands)
            self.y_repetitions = cube_rows
            self.x_repetitions = cube_cols
            # Expand (base_bands) -> (base_bands,1,1)
            data = np.expand_dims(base, axis=(1,2))

        elif base_dimensions == 2:
            # 2D: Spatial Input
            (base_rows, base_cols) = np.shape(base)
            self.t_repetitions = cube_bands
            self.y_repetitions = divide_evenly(cube_rows, base_rows)
            self.x_repetitions = divide_evenly(cube_cols, base_cols)
            # Expand (base_rows,base_cols) -> (1,base_rows,base_cols)
            data = np.expand_dims(base, axis=0)

        elif base_dimensions == 3:
            # 3D: Spatio-Temporal Input
            (base_bands, base_rows, base_cols) = np.shape(base)
            self.t_repetitions = divide_evenly(cube_bands, base_bands)
            self.y_repetitions = divide_evenly(cube_rows, base_rows)
            self.x_repetitions = divide_evenly(cube_cols, base_cols)
            data = np.asarray(base)

        else:
            # 4D+: Invalid Input
            raise ValueError("Invalid input: base must have 0-3 dimensions.")
        self.data = data

    def get(self, t, y, x):
        """
        Return the scalar value at index (t,y,x) by translating these cube coordinates
        to base coordinates and looking up the value within the base data.

        NOTE: Indices may be negative.
        """
        # Select value by spatio-temporal coordinate
        return self.data[t // self.t_repetitions,
                         y // self.y_repetitions,
                         x // self.x_repetitions]

    @cy.profile(False)
    @cy.cdivision(True)
    #@cy.nogil
    @cy.exceptval(check=False)
    @cy.boundscheck(False)
    def get_float(self, t, y, x):
        return self.data[
            t // self.t_repetitions,
            y // self.y_repetitions,
            x // self.x_repetitions
            ]


    def getTimeSeries(self, t_range, y, x):
        """
        Return the 1D array given by the slice (t_range,y,x) by translating these cube
        coordinates to base coordinates, looking up the array slice within the base data,
        and expanding it back to the cube_shape resolution.

        NOTE: Indices may be negative.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument range
        (t_start, t_stop_exclusive) = to_positive_index_range(t_range, self.shape[0])
        t_stop = t_stop_exclusive - 1
        # Translate high-res coordinates to low-res coordinates
        t_start_chunk = t_start // self.t_repetitions
        t_stop_chunk  = t_stop  // self.t_repetitions
        y_chunk       = y       // self.y_repetitions
        x_chunk       = x       // self.x_repetitions
        # Select the array slice that completely contains all low-res coordinates
        low_res_time = self.data[t_start_chunk:(t_stop_chunk + 1),
                                 y_chunk,
                                 x_chunk]
        # Expand the low-res slice into a high-res slice
        high_res_time = maybe_repeat_array(low_res_time, (0, self.t_repetitions))
        # Translate high-res global coordinates to high-res slice coordinates
        t_chunk_origin = t_start_chunk * self.t_repetitions
        t_start_idx    = t_start - t_chunk_origin
        t_stop_idx     = t_stop  - t_chunk_origin
        # Select the array slice that matches the high-res slice coordinates
        return high_res_time[t_start_idx:(t_stop_idx + 1)]


    def getSpatialPlane(self, t, y_range, x_range):
        """
        Return the 2D array given by the slice (t,y_range,x_range) by translating these
        cube coordinates to base coordinates, looking up the array slice within the base
        data, and expanding it back to the cube_shape resolution.

        NOTE: Indices may be negative.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument ranges
        (y_start, y_stop_exclusive) = to_positive_index_range(y_range, self.shape[1])
        (x_start, x_stop_exclusive) = to_positive_index_range(x_range, self.shape[2])
        y_stop = y_stop_exclusive - 1
        x_stop = x_stop_exclusive - 1
        # Translate high-res coordinates to low-res coordinates
        t_chunk       = t       // self.t_repetitions
        y_start_chunk = y_start // self.y_repetitions
        y_stop_chunk  = y_stop  // self.y_repetitions
        x_start_chunk = x_start // self.x_repetitions
        x_stop_chunk  = x_stop  // self.x_repetitions
        # Select the array slice that completely contains all low-res coordinates
        low_res_space = self.data[t_chunk,
                                  y_start_chunk:(y_stop_chunk + 1),
                                  x_start_chunk:(x_stop_chunk + 1)]
        # Expand the low-res slice into a high-res slice
        high_res_space = reduce(maybe_repeat_array,
                                ((0, self.y_repetitions),
                                 (1, self.x_repetitions)),
                                low_res_space)
        # Translate high-res global coordinates to high-res slice coordinates
        y_chunk_origin = y_start_chunk * self.y_repetitions
        x_chunk_origin = x_start_chunk * self.x_repetitions
        y_start_idx    = y_start - y_chunk_origin
        y_stop_idx     = y_stop  - y_chunk_origin
        x_start_idx    = x_start - x_chunk_origin
        x_stop_idx     = x_stop  - x_chunk_origin
        # Select the array slice that matches the high-res slice coordinates
        return high_res_space[y_start_idx:(y_stop_idx + 1),
                              x_start_idx:(x_stop_idx + 1)]


    def getSubcube(self, t_range, y_range, x_range):
        """
        Return the 3D array given by the slice (t_range,y_range,x_range) by translating
        these cube coordinates to base coordinates, looking up the array slice within the
        base data, and expanding it back to the cube_shape resolution.

        NOTE: Indices may be negative.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument ranges
        (t_start, t_stop_exclusive) = to_positive_index_range(t_range, self.shape[0])
        (y_start, y_stop_exclusive) = to_positive_index_range(y_range, self.shape[1])
        (x_start, x_stop_exclusive) = to_positive_index_range(x_range, self.shape[2])
        t_stop = t_stop_exclusive - 1
        y_stop = y_stop_exclusive - 1
        x_stop = x_stop_exclusive - 1
        # Translate high-res coordinates to low-res coordinates
        t_start_chunk = t_start // self.t_repetitions
        t_stop_chunk  = t_stop  // self.t_repetitions
        y_start_chunk = y_start // self.y_repetitions
        y_stop_chunk  = y_stop  // self.y_repetitions
        x_start_chunk = x_start // self.x_repetitions
        x_stop_chunk  = x_stop  // self.x_repetitions
        # Select the array slice that completely contains all low-res coordinates
        low_res_cube = self.data[t_start_chunk:(t_stop_chunk + 1),
                                 y_start_chunk:(y_stop_chunk + 1),
                                 x_start_chunk:(x_stop_chunk + 1)]
        # Expand the low-res slice into a high-res slice
        high_res_cube = reduce(maybe_repeat_array,
                               ((0, self.t_repetitions),
                                (1, self.y_repetitions),
                                (2, self.x_repetitions)),
                               low_res_cube)
        # Translate high-res global coordinates to high-res slice coordinates
        t_chunk_origin = t_start_chunk * self.t_repetitions
        y_chunk_origin = y_start_chunk * self.y_repetitions
        x_chunk_origin = x_start_chunk * self.x_repetitions
        t_start_idx    = t_start - t_chunk_origin
        t_stop_idx     = t_stop  - t_chunk_origin
        y_start_idx    = y_start - y_chunk_origin
        y_stop_idx     = y_stop  - y_chunk_origin
        x_start_idx    = x_start - x_chunk_origin
        x_stop_idx     = x_stop  - x_chunk_origin
        # Select the array slice that matches the high-res slice coordinates
        return high_res_cube[t_start_idx:(t_stop_idx + 1),
                             y_start_idx:(y_stop_idx + 1),
                             x_start_idx:(x_stop_idx + 1)]


    def __getFullyRealizedCube(self):
        """
        Return the 3D array created by expanding the base data to the cube_shape resolution.
        Wherever possible, Numpy broadcasting is used to avoid memory allocation along
        constant array dimensions.
        """
        base_dimensions = np.ndim(self.base)

        if base_dimensions == 0:
            # 0D: Constant Input
            # Broadcast (0,0,0) -> (t,y,x)
            return np.broadcast_to(self.data, self.shape)

        elif base_dimensions == 1:
            # 1D: Time-Series Input
            # Repeat (t0,1,1) -> (t,1,1)
            repeated_array = maybe_repeat_array(self.data, (0, self.t_repetitions))
            # Broadcast (t,1,1) -> (t,y,x)
            return np.broadcast_to(repeated_array, self.shape)

        elif base_dimensions == 2:
            # 2D: Spatial Input
            # Repeat (1,y0,x0) -> (1,y,x)
            repeated_array = reduce(maybe_repeat_array,
                                    ((1, self.y_repetitions),
                                     (2, self.x_repetitions)),
                                    self.data)
            # Broadcast (1,y,x) -> (t,y,x)
            return np.broadcast_to(repeated_array, self.shape)

        else:
            # 3D: Spatio-Temporal Input
            # Repeat (t0,y0,x0) -> (t,y,x)
            return reduce(maybe_repeat_array,
                          ((0, self.t_repetitions),
                           (1, self.y_repetitions),
                           (2, self.x_repetitions)),
                          self.data)


    def getFullyRealizedCube(self, cache=False):
        """
        Return the 3D array created by expanding the base data to the cube_shape resolution.
        Wherever possible, Numpy broadcasting is used to avoid memory allocation along
        constant array dimensions. When cache == True, this expanded 3D array is cached
        within the SpaceTimeCube object for future immediate retrieval.
        """
        if hasattr(self, "cube"):
            return self.cube
        else:
            cube = self.__getFullyRealizedCube()
            if cache is True:
                self.cube = cube
            return cube


    def releaseFullyRealizedCube(self):
        """
        Deletes the cached fully realized cube if it exists.
        """
        if hasattr(self, "cube"):
            delattr(self, "cube")
# space-time-cube-class ends here
# [[file:../../org/pyretechnics.org::lazy-space-time-cube-class][lazy-space-time-cube-class]]
class LazySpaceTimeCube:
    """
    Create an object that represents a 3D array with dimensions (T,Y,X) given by cube_shape.
    Internally, data is stored as an initially empty 3D array of SpaceTimeCube objects.
    Whenever a point value or contiguous space-time region of values is requested, identify
    which SpaceTimeCubes contain the requested coordinates, load them into the cache array
    by calling load_subcube for any that are not already present, request the values from
    these SpaceTimeCubes, combine them together if necessary, and return the resulting scalar
    value or array to the caller.
    """
    def __init__(self, cube_shape, subcube_shape, load_subcube):
        """
        NOTE: The resolutions in cube_shape must be exact multiples of those in subcube_shape.
        """
        # Ensure that cube_shape and subcube_shape both contain 3 values or throw an error
        (cube_bands, cube_rows, cube_cols) = cube_shape
        (subcube_bands, subcube_rows, subcube_cols) = subcube_shape

        # Ensure that cube_shape and subcube_shape only contain positive integers or throw an error
        if not(all(map(is_pos_int, cube_shape + subcube_shape))):
            raise ValueError("The cube_shape and subcube_shape must only contain positive integers.")

        # Ensure that cube_shape is divided evenly by subcube_shape or throw an error
        cache_bands = divide_evenly(cube_bands, subcube_bands)
        cache_rows  = divide_evenly(cube_rows, subcube_rows)
        cache_cols  = divide_evenly(cube_cols, subcube_cols)

        # Store the cube metadata, subcube_shape, cache_shape, cache, and load_subcube functions for later
        self.ndim          = 3
        self.size          = cube_bands * cube_rows * cube_cols
        self.shape         = cube_shape
        self.subcube_shape = subcube_shape
        self.cache_shape   = (cache_bands, cache_rows, cache_cols)
        self.cache         = np.empty(self.cache_shape, dtype=object)
        self.load_subcube  = load_subcube


    def __getOrLoadSubcube(self, cache_t, cache_y, cache_x):
        """
        Return the SpaceTimeCube stored at self.cache[cache_t, cache_y, cache_x] if it
        has already been loaded. Otherwise, call self.load_subcube to load it, store
        it in self.cache, and return it.
        """
        subcube = self.cache[cache_t, cache_y, cache_x]
        if subcube:
            return subcube
        else:
            subcube = self.load_subcube((cache_t, cache_y, cache_x), self.subcube_shape)
            self.cache[cache_t, cache_y, cache_x] = subcube
            return subcube


    def get(self, t, y, x):
        """
        Return the scalar value at index (t,y,x) by translating these cube coordinates
        to cache and subcube coordinates, loading the matching subcube into the cache grid
        if not already present, and looking up the value within this subcube.

        NOTE: Indices may be negative provided that your load_subcube function can handle
              negative indices in its cache_index argument.
        """
        (subcube_bands, subcube_rows, subcube_cols) = self.subcube_shape
        (cache_t, subcube_t) = divmod(t, subcube_bands)
        (cache_y, subcube_y) = divmod(y, subcube_rows)
        (cache_x, subcube_x) = divmod(x, subcube_cols)
        subcube = self.__getOrLoadSubcube(cache_t, cache_y, cache_x)
        return subcube.get(subcube_t, subcube_y, subcube_x)


    def getTimeSeries(self, t_range, y, x):
        """
        Return the 1D array given by the slice (t_range,y,x) by translating these cube
        coordinates to cache and subcube coordinates, loading the matching subcubes into
        the cache grid if not already present, looking up the array slices within each
        subcube, and merging them together into a single 1D array.

        NOTE: Indices may be negative provided that your load_subcube function can handle
              negative indices in its cache_index argument.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument range
        (t_start, t_stop_exclusive) = to_positive_index_range(t_range, self.shape[0])
        t_stop = t_stop_exclusive - 1
        # Translate high-res coordinates to cache and subcube coordinates
        (subcube_bands, subcube_rows, subcube_cols) = self.subcube_shape
        (cache_t_start, subcube_t_start) = divmod(t_start, subcube_bands)
        (cache_t_stop,  subcube_t_stop)  = divmod(t_stop,  subcube_bands)
        (cache_y,       subcube_y)       = divmod(y,       subcube_rows)
        (cache_x,       subcube_x)       = divmod(x,       subcube_cols)
        # Load, expand, and combine subcubes
        return np.concatenate(
            [self.__getOrLoadSubcube(cache_t,
                                     cache_y,
                                     cache_x
                                    ).getTimeSeries(
                                        (subcube_t_start    if cache_t == cache_t_start else 0,
                                         subcube_t_stop + 1 if cache_t == cache_t_stop  else subcube_bands),
                                        subcube_y,
                                        subcube_x
                                    )
             for cache_t in range(cache_t_start, cache_t_stop + 1)]
        )


    def getSpatialPlane(self, t, y_range, x_range):
        """
        Return the 2D array given by the slice (t,y_range,x_range) by translating these
        cube coordinates to cache and subcube coordinates, loading the matching subcubes
        into the cache grid if not already present, looking up the array slices within each
        subcube, and merging them together into a single 2D array.

        NOTE: Indices may be negative provided that your load_subcube function can handle
              negative indices in its cache_index argument.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument ranges
        (y_start, y_stop_exclusive) = to_positive_index_range(y_range, self.shape[1])
        (x_start, x_stop_exclusive) = to_positive_index_range(x_range, self.shape[2])
        y_stop = y_stop_exclusive - 1
        x_stop = x_stop_exclusive - 1
        # Translate high-res coordinates to cache and subcube coordinates
        (subcube_bands, subcube_rows, subcube_cols) = self.subcube_shape
        (cache_t,       subcube_t)       = divmod(t,       subcube_bands)
        (cache_y_start, subcube_y_start) = divmod(y_start, subcube_rows)
        (cache_y_stop,  subcube_y_stop)  = divmod(y_stop,  subcube_rows)
        (cache_x_start, subcube_x_start) = divmod(x_start, subcube_cols)
        (cache_x_stop,  subcube_x_stop)  = divmod(x_stop,  subcube_cols)
        # Load, expand, and combine subcubes
        return np.block(
            [[self.__getOrLoadSubcube(cache_t,
                                      cache_y,
                                      cache_x
                                      ).getSpatialPlane(
                                          subcube_t,
                                          (subcube_y_start    if cache_y == cache_y_start else 0,
                                           subcube_y_stop + 1 if cache_y == cache_y_stop  else subcube_rows),
                                          (subcube_x_start    if cache_x == cache_x_start else 0,
                                           subcube_x_stop + 1 if cache_x == cache_x_stop  else subcube_cols)
                                      )
              for cache_x in range(cache_x_start, cache_x_stop + 1)]
             for cache_y in range(cache_y_start, cache_y_stop + 1)]
        )


    def getSubcube(self, t_range, y_range, x_range):
        """
        Return the 3D array given by the slice (t_range,y_range,x_range) by translating
        these cube coordinates to cache and subcube coordinates, loading the matching
        subcubes into the cache grid if not already present, looking up the array slices
        within each subcube, and merging them together into a single 3D array.

        NOTE: Indices may be negative provided that your load_subcube function can handle
              negative indices in its cache_index argument.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument ranges
        (t_start, t_stop_exclusive) = to_positive_index_range(t_range, self.shape[0])
        (y_start, y_stop_exclusive) = to_positive_index_range(y_range, self.shape[1])
        (x_start, x_stop_exclusive) = to_positive_index_range(x_range, self.shape[2])
        t_stop = t_stop_exclusive - 1
        y_stop = y_stop_exclusive - 1
        x_stop = x_stop_exclusive - 1
        # Translate high-res coordinates to cache and subcube coordinates
        (subcube_bands, subcube_rows, subcube_cols) = self.subcube_shape
        (cache_t_start, subcube_t_start) = divmod(t_start, subcube_bands)
        (cache_t_stop,  subcube_t_stop)  = divmod(t_stop,  subcube_bands)
        (cache_y_start, subcube_y_start) = divmod(y_start, subcube_rows)
        (cache_y_stop,  subcube_y_stop)  = divmod(y_stop,  subcube_rows)
        (cache_x_start, subcube_x_start) = divmod(x_start, subcube_cols)
        (cache_x_stop,  subcube_x_stop)  = divmod(x_stop,  subcube_cols)
        # Load, expand, and combine subcubes
        return np.block(
            [[[self.__getOrLoadSubcube(cache_t,
                                       cache_y,
                                       cache_x
                                       ).getSubcube(
                                           (subcube_t_start    if cache_t == cache_t_start else 0,
                                            subcube_t_stop + 1 if cache_t == cache_t_stop  else subcube_bands),
                                           (subcube_y_start    if cache_y == cache_y_start else 0,
                                            subcube_y_stop + 1 if cache_y == cache_y_stop  else subcube_rows),
                                           (subcube_x_start    if cache_x == cache_x_start else 0,
                                            subcube_x_stop + 1 if cache_x == cache_x_stop  else subcube_cols)
                                       )
               for cache_x in range(cache_x_start, cache_x_stop + 1)]
              for cache_y in range(cache_y_start, cache_y_stop + 1)]
             for cache_t in range(cache_t_start, cache_t_stop + 1)]
        )


    def getFullyRealizedCube(self, cache=False):
        raise ValueError("getFullyRealizedCube is not implemented for LazySpaceTimeCube.\n"
                         + "You probably don't want to do this anyway.")


    def releaseFullyRealizedCube(self):
        raise ValueError("releaseFullyRealizedCube is not implemented for LazySpaceTimeCube.\n"
                         + "You probably don't want to do this anyway.")
# lazy-space-time-cube-class ends here
