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


def maybe_repeat_array(array, axis_repetitions):
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
class SpaceTimeCube:
    """
    TODO: Add docstring.
    """
    def __init__(self, cube_shape, base):
        """
        NOTE: (t,y,x) = (0,0,0) is the upper-left corner of the cube in the first timestep.
        NOTE: cube_shape >= base
        """
        # Ensure that cube_shape contains 3 values or throw an error
        (cube_bands, cube_rows, cube_cols) = cube_shape

        # Ensure that cube_shape only contains positive integers or throw an error
        if not(all(map(is_pos_int, cube_shape))):
            raise ValueError("The cube_shape must only contain positive integers.")

        # Store the cube metadata for later
        self.ndim  = 3
        self.size  = cube_bands * cube_rows * cube_cols
        self.shape = cube_shape
        self.base  = base

        match np.ndim(base):
            # 0D: Constant Input
            case 0:
                self.t_repetitions = cube_bands
                self.y_repetitions = cube_rows
                self.x_repetitions = cube_cols
                self.data = np.asarray([[[base]]])

            # 1D: Time-Series Input
            case 1:
                base_bands = len(base)
                self.t_repetitions = divide_evenly(cube_bands, base_bands)
                self.y_repetitions = cube_rows
                self.x_repetitions = cube_cols
                # Expand (base_bands) -> (base_bands,1,1)
                self.data = np.expand_dims(base, axis=(1,2))

            # 2D: Spatial Input
            case 2:
                (base_rows, base_cols) = np.shape(base)
                self.t_repetitions = cube_bands
                self.y_repetitions = divide_evenly(cube_rows, base_rows)
                self.x_repetitions = divide_evenly(cube_cols, base_cols)
                # Expand (base_rows,base_cols) -> (1,base_rows,base_cols)
                self.data = np.expand_dims(base, axis=0)

            # 3D: Spatio-Temporal Input
            case 3:
                (base_bands, base_rows, base_cols) = np.shape(base)
                self.t_repetitions = divide_evenly(cube_bands, base_bands)
                self.y_repetitions = divide_evenly(cube_rows, base_rows)
                self.x_repetitions = divide_evenly(cube_cols, base_cols)
                self.data = np.asarray(base)

            # 4D+: Invalid Input
            case _:
                raise ValueError("Invalid input: base must have 0-3 dimensions.")


    def get(self, t, y, x):
        """
        Returns a scalar value.
        """
        # Select value by spatio-temporal coordinate
        return self.data[t // self.t_repetitions,
                         y // self.y_repetitions,
                         x // self.x_repetitions]


    # NOTE: Ranges provide inclusion:inclusion semantics
    # NOTE: None and -1 cannot be passed in
    def getTimeSeries(self, t_range, y, x):
        """
        Returns a 1D array.
        """
        # Destructure the argument range
        (t_start, t_stop) = t_range
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


    # NOTE: Ranges provide inclusion:inclusion semantics
    # NOTE: None and -1 cannot be passed in
    def getSpatialPlane(self, t, y_range, x_range):
        """
        Returns a 2D array.
        """
        # Destructure the argument ranges
        (y_start, y_stop) = y_range
        (x_start, x_stop) = x_range
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


    # NOTE: Ranges provide inclusion:inclusion semantics
    # NOTE: None and -1 cannot be passed in
    def getSubcube(self, t_range, y_range, x_range):
        """
        Returns a 3D array.
        """
        # Destructure the argument ranges
        (t_start, t_stop) = t_range
        (y_start, y_stop) = y_range
        (x_start, x_stop) = x_range
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
        match np.ndim(self.base):
            # 0D: Constant Input
            case 0:
                # Broadcast (0,0,0) -> (t,y,x)
                return np.broadcast_to(self.data, self.shape)

            # 1D: Time-Series Input
            case 1:
                # Repeat (t0,1,1) -> (t,1,1)
                repeated_array = maybe_repeat_array(self.data, (0, self.t_repetitions))
                # Broadcast (t,1,1) -> (t,y,x)
                return np.broadcast_to(repeated_array, self.shape)

            # 2D: Spatial Input
            case 2:
                # Repeat (1,y0,x0) -> (1,y,x)
                repeated_array = reduce(maybe_repeat_array,
                                        ((1, self.y_repetitions),
                                         (2, self.x_repetitions)),
                                        self.data)
                # Broadcast (1,y,x) -> (t,y,x)
                return np.broadcast_to(repeated_array, self.shape)

            # 3D: Spatio-Temporal Input
            case 3:
                # Repeat (t0,y0,x0) -> (t,y,x)
                return reduce(maybe_repeat_array,
                              ((0, self.t_repetitions),
                               (1, self.y_repetitions),
                               (2, self.x_repetitions)),
                              self.data)


    def getFullyRealizedCube(self, cache=False):
        """
        Returns a 3D array.
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
    TODO: Add docstring.
    """
    def __init__(self, cube_shape, subcube_shape, load_subcube):
        """
        NOTE: (t,y,x) = (0,0,0) is the upper-left corner of the cube in the first timestep.
        NOTE: cube_shape >= subcube_shape
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


    def __getOrLoadSubcube(cache_t, cache_y, cache_x):
        """
        Returns a SpaceTimeCube.
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
        Returns a scalar value.
        """
        (subcube_bands, subcube_rows, subcube_cols) = self.subcube_shape
        (cache_t, subcube_t) = divmod(t, subcube_bands)
        (cache_y, subcube_y) = divmod(y, subcube_rows)
        (cache_x, subcube_x) = divmod(x, subcube_cols)
        subcube = self.__getOrLoadSubcube(cache_t, cache_y, cache_x)
        return subcube.get(subcube_t, subcube_y, subcube_x)


    # NOTE: Ranges provide inclusion:inclusion semantics
    # NOTE: None and -1 cannot be passed in
    def getTimeSeries(self, t_range, y, x):
        """
        Returns a 1D array.
        """
        # Destructure the argument range
        (t_start, t_stop) = t_range
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
                                        (subcube_t_start if cache_t == cache_t_start else 0,
                                         subcube_t_stop  if cache_t == cache_t_stop  else subcube_bands - 1),
                                        subcube_y,
                                        subcube_x
                                    )
             for cache_t in range(cache_t_start, cache_t_stop + 1)]
        )


    # NOTE: Ranges provide inclusion:inclusion semantics
    # NOTE: None and -1 cannot be passed in
    def getSpatialPlane(self, t, y_range, x_range):
        """
        Returns a 2D array.
        """
        # Destructure the argument ranges
        (y_start, y_stop) = y_range
        (x_start, x_stop) = x_range
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
                                          (subcube_y_start if cache_y == cache_y_start else 0,
                                           subcube_y_stop  if cache_y == cache_y_stop  else subcube_rows - 1),
                                          (subcube_x_start if cache_x == cache_x_start else 0,
                                           subcube_x_stop  if cache_x == cache_x_stop  else subcube_cols - 1)
                                      )
              for cache_x in range(cache_x_start, cache_x_stop + 1)]
             for cache_y in range(cache_y_start, cache_y_stop + 1)]
        )


    # NOTE: Ranges provide inclusion:inclusion semantics
    # NOTE: None and -1 cannot be passed in
    def getSubcube(self, t_range, y_range, x_range):
        """
        Returns a 3D array.
        """
        # Destructure the argument ranges
        (t_start, t_stop) = t_range
        (y_start, y_stop) = y_range
        (x_start, x_stop) = x_range
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
                                           (subcube_t_start if cache_t == cache_t_start else 0,
                                            subcube_t_stop  if cache_t == cache_t_stop  else subcube_bands - 1),
                                           (subcube_y_start if cache_y == cache_y_start else 0,
                                            subcube_y_stop  if cache_y == cache_y_stop  else subcube_rows - 1),
                                           (subcube_x_start if cache_x == cache_x_start else 0,
                                            subcube_x_stop  if cache_x == cache_x_stop  else subcube_cols - 1)
                                       )
               for cache_x in range(cache_x_start, cache_x_stop + 1)]
              for cache_y in range(cache_y_start, cache_y_stop + 1)]
             for cache_t in range(cache_t_start, cache_t_stop + 1)]
        )


    def getFullyRealizedCube(self, cache=False):
        raise ValueError("getFullyRealizedCube is not implemented for LazySpaceTimeCube.\n"
                         + "You probably don't want to do this anyway.")
# lazy-space-time-cube-class ends here
