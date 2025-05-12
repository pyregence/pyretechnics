# [[file:../../org/pyretechnics.org::space-time-cube-imports][space-time-cube-imports]]
import cython
import cython as cy
from functools import reduce
import numpy as np
if cython.compiled:
    from cython.cimports.numpy import ndarray
    from cython.cimports.pyretechnics.cy_types import pyidx
else:
    from numpy import ndarray
    from pyretechnics.py_types import pyidx
# space-time-cube-imports ends here
# [[file:../../org/pyretechnics.org::space-time-cube-utilities][space-time-cube-utilities]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def is_pos_int(x: object) -> cy.bint:
    return isinstance(x, int) and x > 0


@cy.cfunc
def divide_evenly(dividend: cy.int, divisor: cy.int) -> cy.int:
    if divisor == 0:
        raise ValueError(f"{divisor} may not be zero.")
    else:
        quotient : cy.int = dividend // divisor
        remainder: cy.int = dividend % divisor
        if remainder == 0:
            return quotient
        else:
            raise ValueError(f"{dividend} must be an exact multiple of {divisor}.")


@cy.ccall
@cy.exceptval(check=False)
def to_positive_index_range(index_range: tuple[pyidx, pyidx]|None, axis_length: pyidx) -> tuple[pyidx, pyidx]:
    """
    Translate None and negative indices to positive indices.
    """
    if index_range is None:
        return (0, axis_length)
    else:
        start: pyidx|None = index_range[0]
        stop : pyidx|None = index_range[1]
        return (
            0 if start is None else axis_length + start if start < 0 else start,
            axis_length if stop is None else axis_length + stop if stop < 0 else stop
        )


@cy.cfunc
def stretch_array(old_array: ndarray, new_length: cy.int, repetitions: cy.float) -> ndarray:
    new_array: ndarray = np.zeros(new_length, dtype=old_array.dtype)
    i        : pyidx
    for i in range(new_length):
        new_array[i] = old_array[int(i / repetitions)]
    return new_array


@cy.cfunc
def maybe_repeat_array(maybe_array: ndarray, axis_repetitions: tuple[pyidx, cy.float]) -> ndarray:
    """
    Return a new array that is created by repeating the elements from the input
    array repetitions times along the specified array axis. Avoid allocating
    new memory if repetitions == 1 or if the repeated array axis has length 1.
    """
    axis       : pyidx    = axis_repetitions[0]
    repetitions: cy.float = axis_repetitions[1]
    old_array  : ndarray  = np.asarray(maybe_array)
    array_dims : pyidx    = old_array.ndim
    array_shape: list     = list(np.shape(old_array))
    axis_length: cy.int   = array_shape[axis]
    if repetitions == 1.0:
        # no repetitions necessary
        return old_array
    elif axis_length == 1:
        # broadcast single-element axis repetitions times
        array_shape[axis] = int(repetitions)
        return np.broadcast_to(old_array, array_shape)
    elif repetitions % 1.0 == 0.0:
        # repeat each element on the chosen axis repetitions times
        return np.repeat(old_array, int(repetitions), axis)
    elif axis == 0 and array_dims == 1:
        # populate a new 1D array of the expected length by translating its indices into the original array
        new_axis_length: cy.int = int(axis_length * repetitions)
        return stretch_array(old_array, new_axis_length, repetitions)
    elif axis == 0 and array_dims == 2:
        # populate a new 2D array of the expected length by translating its indices into the original array
        new_rows: cy.int = int(axis_length * repetitions)
        return np.stack([old_array[int(i / repetitions)] for i in range(new_rows)])
    elif axis == 1 and array_dims == 2:
        # populate a new 2D array of the expected length by translating its indices into the original array
        old_rows: cy.int = array_shape[0]
        new_cols: cy.int = int(axis_length * repetitions)
        return np.stack([stretch_array(old_array[i], new_cols, repetitions) for i in range(old_rows)])
    elif axis == 0 and array_dims == 3:
        # populate a new 3D array of the expected length by translating its indices into the original array
        new_bands: cy.int = int(axis_length * repetitions)
        return np.stack([old_array[int(b / repetitions)] for b in range(new_bands)])
    elif axis == 1 and array_dims == 3:
        # populate a new 2D array of the expected length by translating its indices into the original array
        old_bands: cy.int = array_shape[0]
        new_rows : cy.int = int(axis_length * repetitions)
        return np.stack([
            np.stack([
                old_array[b, int(i / repetitions)]
                for i in range(new_rows)
            ])
            for b in range(old_bands)
        ])
    elif axis == 2 and array_dims == 3:
        # populate a new 3D array of the expected length by translating its indices into the original array
        old_bands: cy.int = array_shape[0]
        old_rows : cy.int = array_shape[1]
        new_cols : cy.int = int(axis_length * repetitions)
        return np.stack([
            np.stack([
                stretch_array(old_array[b,i], new_cols, repetitions)
                for i in range(old_rows)
            ])
            for b in range(old_bands)
        ])
    else:
        raise ValueError("Floating point repetitions are only supported for 1D, 2D, and 3D arrays.")
# space-time-cube-utilities ends here
# [[file:../../org/pyretechnics.org::ispace-time-cube-class][ispace-time-cube-class]]
@cy.cclass
class ISpaceTimeCube:
    @cy.cfunc
    @cy.exceptval(check=False)
    def get(self, t: pyidx, y: pyidx, x: pyidx) -> cy.float:
        pass
# ispace-time-cube-class ends here
# [[file:../../org/pyretechnics.org::space-time-cube-class][space-time-cube-class]]
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
    ndim         : cy.int
    size         : cy.ulonglong
    shape        : tuple[cy.int, cy.int, cy.int]
    base         : object
    t_repetitions: cy.float
    y_repetitions: cy.float
    x_repetitions: cy.float
    data         : cy.float[:,:,::1] # FIXME: Restore polymorphism for the underlying Numpy arrays


    def __init__(self, cube_shape: tuple[int, int, int], base: object) -> cy.void:
        """
        NOTE: The resolutions in cube_shape must be exact multiples of any existing dimensions
              in the base data. If base is not a Numpy float32 array, a new array will be allocated.
        """
        # Ensure that cube_shape contains 3 values or throw an error
        if len(cube_shape) != 3:
            raise ValueError("The cube_shape must contain exactly three values.")

        # Unpack the cube_shape values without type-checking
        cube_bands_: object = cube_shape[0]
        cube_rows_ : object = cube_shape[1]
        cube_cols_ : object = cube_shape[2]

        # Ensure that cube_shape only contains positive integers or throw an error
        if not(is_pos_int(cube_bands_) and is_pos_int(cube_rows_) and is_pos_int(cube_cols_)):
            raise ValueError("The cube_shape must contain only positive integers.")

        # Cast the cube_shape values as primitive ints
        cube_bands: cy.int = cube_bands_
        cube_rows : cy.int = cube_rows_
        cube_cols : cy.int = cube_cols_

        # Store the cube metadata for later
        self.ndim  = 3
        self.size  = cube_bands * cube_rows * cube_cols
        self.shape = (cube_bands, cube_rows, cube_cols)
        self.base  = base

        # Store the base data as a 3D array along with its axis repetitions
        base_dimensions: cy.int = np.ndim(base)

        if base_dimensions == 0:
            # 0D: Constant Input
            self.t_repetitions = cube_bands
            self.y_repetitions = cube_rows
            self.x_repetitions = cube_cols
            self.data          = np.asarray([[[base]]], dtype=np.float32)

        elif base_dimensions == 1:
            # 1D: Time-Series Input
            base_bands: cy.int = len(base)
            self.t_repetitions = float(cube_bands) / float(base_bands)
            self.y_repetitions = cube_rows
            self.x_repetitions = cube_cols
            # Ensure that the cube_shape is not smaller than the base shape
            if cube_bands < base_bands:
                raise ValueError("The cube_shape may not be smaller than the base shape.")
            # Warn if any repetitions are not whole numbers
            if self.t_repetitions % 1.0 != 0.0:
                print("WARNING: Input data's shape does not evenly divide the cube_shape."
                      + " Index lookups beyond the cube's edge may mistakenly return values without errors.")
            # Warn if base is not a Numpy float32 array
            if not(isinstance(base, np.ndarray)) or (base.dtype != np.float32):
                print("WARNING: Input data is not a Numpy float32 array. Data will be copied into SpaceTimeCube.",
                      flush=True)
            # Expand (base_bands) -> (base_bands,1,1)
            self.data = np.expand_dims(np.asarray(base, dtype=np.float32), axis=(1,2))

        elif base_dimensions == 2:
            # 2D: Spatial Input
            base_shape: tuple  = np.shape(base)
            base_rows : cy.int = base_shape[0]
            base_cols : cy.int = base_shape[1]
            self.t_repetitions = cube_bands
            self.y_repetitions = float(cube_rows) / float(base_rows)
            self.x_repetitions = float(cube_cols) / float(base_cols)
            # Ensure that the cube_shape is not smaller than the base shape
            if cube_rows < base_rows or cube_cols < base_cols:
                raise ValueError("The cube_shape may not be smaller than the base shape.")
            # Warn if any repetitions are not whole numbers
            if self.y_repetitions % 1.0 != 0.0 or self.x_repetitions % 1.0 != 0.0:
                print("WARNING: Input data's shape does not evenly divide the cube_shape."
                      + " Index lookups beyond the cube's edge may mistakenly return values without errors.")
            # Warn if base is not a Numpy float32 array
            if not(isinstance(base, np.ndarray)) or (base.dtype != np.float32):
                print("WARNING: Input data is not a Numpy float32 array. Data will be copied into SpaceTimeCube.",
                      flush=True)
            # Expand (base_rows,base_cols) -> (1,base_rows,base_cols)
            self.data = np.expand_dims(np.asarray(base, dtype=np.float32), axis=0)

        elif base_dimensions == 3:
            # 3D: Spatio-Temporal Input
            base_shape: tuple  = np.shape(base)
            base_bands: cy.int = base_shape[0]
            base_rows : cy.int = base_shape[1]
            base_cols : cy.int = base_shape[2]
            self.t_repetitions = float(cube_bands) / float(base_bands)
            self.y_repetitions = float(cube_rows) / float(base_rows)
            self.x_repetitions = float(cube_cols) / float(base_cols)
            # Ensure that the cube_shape is not smaller than the base shape
            if cube_bands < base_bands or cube_rows < base_rows or cube_cols < base_cols:
                raise ValueError("The cube_shape may not be smaller than the base shape.")
            # Warn if any repetitions are not whole numbers
            if self.t_repetitions % 1.0 != 0.0 or self.y_repetitions % 1.0 != 0.0 or self.x_repetitions % 1.0 != 0.0:
                print("WARNING: Input data's shape does not evenly divide the cube_shape."
                      + " Index lookups beyond the cube's edge may mistakenly return values without errors.")
            # Warn if base is not a Numpy float32 array
            if not(isinstance(base, np.ndarray)) or (base.dtype != np.float32):
                print("WARNING: Input data is not a Numpy float32 array. Data will be copied into SpaceTimeCube.",
                      flush=True)
            self.data = np.asarray(base, dtype=np.float32)

        else:
            # 4D+: Invalid Input
            raise ValueError("Invalid input: base must have 0-3 dimensions.")


    @cy.ccall
    @cy.exceptval(check=False)
    @cy.boundscheck(True)
    def get(self, t: pyidx, y: pyidx, x: pyidx) -> cy.float:
        """
        Return the scalar value at index (t,y,x) by translating these cube coordinates
        to base coordinates and looking up the value within the base data.

        NOTE: Indices may be negative.
        """
        # Select value by spatio-temporal coordinate
        base_t: pyidx = int(t / self.t_repetitions)
        base_y: pyidx = int(y / self.y_repetitions)
        base_x: pyidx = int(x / self.x_repetitions)
        return self.data[base_t, base_y, base_x]


    @cy.ccall
    @cy.boundscheck(True)
    def getTimeSeries(self, t_range: tuple[pyidx, pyidx]|None, y: pyidx, x: pyidx) -> ndarray:
        """
        Return the 1D array given by the slice (t_range,y,x) by translating these cube
        coordinates to base coordinates, looking up the array slice within the base data,
        and expanding it back to the cube_shape resolution.

        NOTE: Indices may be negative.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument range
        t_range_updated : tuple[pyidx, pyidx] = to_positive_index_range(t_range, self.shape[0])
        t_start         : pyidx               = t_range_updated[0]
        t_stop_exclusive: pyidx               = t_range_updated[1]
        t_stop          : pyidx               = t_stop_exclusive - 1
        # Translate high-res coordinates to low-res coordinates
        t_start_chunk: pyidx = int(t_start / self.t_repetitions)
        t_stop_chunk : pyidx = int(t_stop  / self.t_repetitions)
        y_chunk      : pyidx = int(y       / self.y_repetitions)
        x_chunk      : pyidx = int(x       / self.x_repetitions)
        # Select the array slice that completely contains all low-res coordinates
        low_res_time: ndarray = np.asarray(self.data[t_start_chunk:(t_stop_chunk + 1),
                                                     y_chunk,
                                                     x_chunk])
        # Expand the low-res slice into a high-res slice
        high_res_time: ndarray = maybe_repeat_array(low_res_time, (0, self.t_repetitions))
        # Translate high-res global coordinates to high-res slice coordinates
        t_chunk_origin: pyidx = int(t_start_chunk * self.t_repetitions)
        t_start_idx   : pyidx = t_start - t_chunk_origin
        t_stop_idx    : pyidx = t_stop  - t_chunk_origin
        # Select the array slice that matches the high-res slice coordinates
        return high_res_time[t_start_idx:(t_stop_idx + 1)]


    @cy.ccall
    @cy.boundscheck(True)
    def getSpatialPlane(self,
                        t      : pyidx,
                        y_range: tuple[pyidx, pyidx]|None,
                        x_range: tuple[pyidx, pyidx]|None) -> ndarray:
        """
        Return the 2D array given by the slice (t,y_range,x_range) by translating these
        cube coordinates to base coordinates, looking up the array slice within the base
        data, and expanding it back to the cube_shape resolution.

        NOTE: Indices may be negative.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument ranges
        y_range_updated : tuple[pyidx, pyidx] = to_positive_index_range(y_range, self.shape[1])
        x_range_updated : tuple[pyidx, pyidx] = to_positive_index_range(x_range, self.shape[2])
        y_start         : pyidx               = y_range_updated[0]
        y_stop_exclusive: pyidx               = y_range_updated[1]
        x_start         : pyidx               = x_range_updated[0]
        x_stop_exclusive: pyidx               = x_range_updated[1]
        y_stop          : pyidx               = y_stop_exclusive - 1
        x_stop          : pyidx               = x_stop_exclusive - 1
        # Translate high-res coordinates to low-res coordinates
        t_chunk      : pyidx = int(t       / self.t_repetitions)
        y_start_chunk: pyidx = int(y_start / self.y_repetitions)
        y_stop_chunk : pyidx = int(y_stop  / self.y_repetitions)
        x_start_chunk: pyidx = int(x_start / self.x_repetitions)
        x_stop_chunk : pyidx = int(x_stop  / self.x_repetitions)
        # Select the array slice that completely contains all low-res coordinates
        low_res_space: ndarray = np.asarray(self.data[t_chunk,
                                                      y_start_chunk:(y_stop_chunk + 1),
                                                      x_start_chunk:(x_stop_chunk + 1)])
        # Expand the low-res slice into a high-res slice
        high_res_space: ndarray = reduce(maybe_repeat_array,
                                         ((0, self.y_repetitions),
                                          (1, self.x_repetitions)),
                                         low_res_space)
        # Translate high-res global coordinates to high-res slice coordinates
        y_chunk_origin: pyidx = int(y_start_chunk * self.y_repetitions)
        x_chunk_origin: pyidx = int(x_start_chunk * self.x_repetitions)
        y_start_idx   : pyidx = y_start - y_chunk_origin
        y_stop_idx    : pyidx = y_stop  - y_chunk_origin
        x_start_idx   : pyidx = x_start - x_chunk_origin
        x_stop_idx    : pyidx = x_stop  - x_chunk_origin
        # Select the array slice that matches the high-res slice coordinates
        return high_res_space[y_start_idx:(y_stop_idx + 1),
                              x_start_idx:(x_stop_idx + 1)]


    @cy.ccall
    @cy.boundscheck(True)
    def getSubcube(self,
                   t_range: tuple[pyidx, pyidx]|None,
                   y_range: tuple[pyidx, pyidx]|None,
                   x_range: tuple[pyidx, pyidx]|None) -> ndarray:
        """
        Return the 3D array given by the slice (t_range,y_range,x_range) by translating
        these cube coordinates to base coordinates, looking up the array slice within the
        base data, and expanding it back to the cube_shape resolution.

        NOTE: Indices may be negative.
        NOTE: Range indices may include one or more None values and
              provide (inclusion, exclusion) semantics like Python array slice notation.
        """
        # Destructure the argument ranges
        t_range_updated : tuple[pyidx, pyidx] = to_positive_index_range(t_range, self.shape[0])
        y_range_updated : tuple[pyidx, pyidx] = to_positive_index_range(y_range, self.shape[1])
        x_range_updated : tuple[pyidx, pyidx] = to_positive_index_range(x_range, self.shape[2])
        t_start         : pyidx               = t_range_updated[0]
        t_stop_exclusive: pyidx               = t_range_updated[1]
        y_start         : pyidx               = y_range_updated[0]
        y_stop_exclusive: pyidx               = y_range_updated[1]
        x_start         : pyidx               = x_range_updated[0]
        x_stop_exclusive: pyidx               = x_range_updated[1]
        t_stop          : pyidx               = t_stop_exclusive - 1
        y_stop          : pyidx               = y_stop_exclusive - 1
        x_stop          : pyidx               = x_stop_exclusive - 1
        # Translate high-res coordinates to low-res coordinates
        t_start_chunk: pyidx = int(t_start / self.t_repetitions)
        t_stop_chunk : pyidx = int(t_stop  / self.t_repetitions)
        y_start_chunk: pyidx = int(y_start / self.y_repetitions)
        y_stop_chunk : pyidx = int(y_stop  / self.y_repetitions)
        x_start_chunk: pyidx = int(x_start / self.x_repetitions)
        x_stop_chunk : pyidx = int(x_stop  / self.x_repetitions)
        # Select the array slice that completely contains all low-res coordinates
        low_res_cube: ndarray = np.asarray(self.data[t_start_chunk:(t_stop_chunk + 1),
                                                     y_start_chunk:(y_stop_chunk + 1),
                                                     x_start_chunk:(x_stop_chunk + 1)])
        # Expand the low-res slice into a high-res slice
        high_res_cube: ndarray = reduce(maybe_repeat_array,
                                        ((0, self.t_repetitions),
                                         (1, self.y_repetitions),
                                         (2, self.x_repetitions)),
                                        low_res_cube)
        # Translate high-res global coordinates to high-res slice coordinates
        t_chunk_origin: pyidx = int(t_start_chunk * self.t_repetitions)
        y_chunk_origin: pyidx = int(y_start_chunk * self.y_repetitions)
        x_chunk_origin: pyidx = int(x_start_chunk * self.x_repetitions)
        t_start_idx   : pyidx = t_start - t_chunk_origin
        t_stop_idx    : pyidx = t_stop  - t_chunk_origin
        y_start_idx   : pyidx = y_start - y_chunk_origin
        y_stop_idx    : pyidx = y_stop  - y_chunk_origin
        x_start_idx   : pyidx = x_start - x_chunk_origin
        x_stop_idx    : pyidx = x_stop  - x_chunk_origin
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
            repeated_array = maybe_repeat_array(np.asarray(self.data), (0, self.t_repetitions))
            # Broadcast (t,1,1) -> (t,y,x)
            return np.broadcast_to(repeated_array, self.shape)

        elif base_dimensions == 2:
            # 2D: Spatial Input
            # Repeat (1,y0,x0) -> (1,y,x)
            repeated_array = reduce(maybe_repeat_array,
                                    ((1, self.y_repetitions),
                                     (2, self.x_repetitions)),
                                    np.asarray(self.data))
            # Broadcast (1,y,x) -> (t,y,x)
            return np.broadcast_to(repeated_array, self.shape)

        else:
            # 3D: Spatio-Temporal Input
            # Repeat (t0,y0,x0) -> (t,y,x)
            return reduce(maybe_repeat_array,
                          ((0, self.t_repetitions),
                           (1, self.y_repetitions),
                           (2, self.x_repetitions)),
                          np.asarray(self.data))


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
@cy.cclass
class LazySpaceTimeCube(ISpaceTimeCube):
    """
    Create an object that represents a 3D array with dimensions (T,Y,X) given by cube_shape.
    Internally, data is stored as an initially empty 3D array of SpaceTimeCube objects.
    Whenever a point value or contiguous space-time region of values is requested, identify
    which SpaceTimeCubes contain the requested coordinates, load them into the cache array
    by calling load_subcube for any that are not already present, request the values from
    these SpaceTimeCubes, combine them together if necessary, and return the resulting scalar
    value or array to the caller.
    """
    ndim         : cy.int
    size         : cy.ulonglong
    shape        : tuple[cy.int, cy.int, cy.int]
    subcube_shape: tuple[cy.int, cy.int, cy.int]
    cache_shape  : tuple[cy.int, cy.int, cy.int]
    cache        : ndarray
    load_subcube : object


    def __init__(self,
                 cube_shape   : tuple[int, int, int],
                 subcube_shape: tuple[int, int, int],
                 load_subcube : object) -> cy.void:
        """
        NOTE: The resolutions in cube_shape must be exact multiples of those in subcube_shape.
        """
        # Ensure that cube_shape and subcube_shape both contain 3 values or throw an error
        if len(cube_shape) != 3:
            raise ValueError("The cube_shape must contain exactly three values.")

        if len(subcube_shape) != 3:
            raise ValueError("The subcube_shape must contain exactly three values.")

        # Unpack the cube_shape values without type-checking
        cube_bands_: object = cube_shape[0]
        cube_rows_ : object = cube_shape[1]
        cube_cols_ : object = cube_shape[2]

        # Unpack the subcube_shape values without type-checking
        subcube_bands_: object = subcube_shape[0]
        subcube_rows_ : object = subcube_shape[1]
        subcube_cols_ : object = subcube_shape[2]

        # Ensure that cube_shape and subcube_shape only contain positive integers or throw an error
        if not(is_pos_int(cube_bands_) and is_pos_int(cube_rows_) and is_pos_int(cube_cols_)):
            raise ValueError("The cube_shape must contain only positive integers.")

        if not(is_pos_int(subcube_bands_) and is_pos_int(subcube_rows_) and is_pos_int(subcube_cols_)):
            raise ValueError("The subcube_shape must contain only positive integers.")

        # Cast the cube_shape values as primitive ints
        cube_bands: cy.int = cube_bands_
        cube_rows : cy.int = cube_rows_
        cube_cols : cy.int = cube_cols_

        # Cast the subcube_shape values as primitive ints
        subcube_bands: cy.int = subcube_bands_
        subcube_rows : cy.int = subcube_rows_
        subcube_cols : cy.int = subcube_cols_

        # Ensure that cube_shape is divided evenly by subcube_shape or throw an error
        cache_bands: cy.int = divide_evenly(cube_bands, subcube_bands)
        cache_rows : cy.int = divide_evenly(cube_rows, subcube_rows)
        cache_cols : cy.int = divide_evenly(cube_cols, subcube_cols)

        # Store the cube metadata, subcube_shape, cache_shape, cache, and load_subcube functions for later
        self.ndim          = 3
        self.size          = cube_bands * cube_rows * cube_cols
        self.shape         = (cube_bands, cube_rows, cube_cols)
        self.subcube_shape = (subcube_bands, subcube_rows, subcube_cols)
        self.cache_shape   = (cache_bands, cache_rows, cache_cols)
        self.cache         = np.empty(self.cache_shape, dtype=object)
        self.load_subcube  = load_subcube


    @cy.cfunc
    def __getOrLoadSubcube(self, cache_t: pyidx, cache_y: pyidx, cache_x: pyidx) -> SpaceTimeCube:
        """
        Return the SpaceTimeCube stored at self.cache[cache_t, cache_y, cache_x] if it
        has already been loaded. Otherwise, call self.load_subcube to load it, store
        it in self.cache, and return it.
        """
        subcube: SpaceTimeCube = cy.cast(SpaceTimeCube, self.cache[cache_t, cache_y, cache_x])
        if subcube:
            return subcube
        else:
            subcube = self.load_subcube((cache_t, cache_y, cache_x), self.subcube_shape)
            self.cache[cache_t, cache_y, cache_x] = subcube
            return subcube


    @cy.ccall
    @cy.exceptval(check=False)
    def get(self, t: pyidx, y: pyidx, x: pyidx) -> cy.float:
        """
        Return the scalar value at index (t,y,x) by translating these cube coordinates
        to cache and subcube coordinates, loading the matching subcube into the cache grid
        if not already present, and looking up the value within this subcube.

        NOTE: Indices may be negative provided that your load_subcube function can handle
              negative indices in its cache_index argument.
        """
        # Grab the subcube_shape tuple
        subcube_shape: tuple[cy.int, cy.int, cy.int] = self.subcube_shape

        # Unpack the subcube_shape values
        subcube_bands: cy.int = subcube_shape[0]
        subcube_rows : cy.int = subcube_shape[1]
        subcube_cols : cy.int = subcube_shape[2]

        # Calculate the cache index
        cache_t: pyidx = t // subcube_bands
        cache_y: pyidx = y // subcube_rows
        cache_x: pyidx = x // subcube_cols

        # Calculate the subcube index
        subcube_t: pyidx = t % subcube_bands
        subcube_y: pyidx = y % subcube_rows
        subcube_x: pyidx = x % subcube_cols

        # Fetch the subcube from the cache
        subcube: SpaceTimeCube = self.__getOrLoadSubcube(cache_t, cache_y, cache_x)

        # Look up the scalar value in the subcube at the subcube index
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
