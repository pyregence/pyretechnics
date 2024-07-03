# [[file:../../org/pyretechnics.org::space-time-cube-class][space-time-cube-class]]
from functools import reduce
import numpy as np

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
    return np.repeat(array, repetitions, axis) if (repetitions > 1) else array


class SpaceTimeCube:
    """
    TODO: Add docstring.
    """
    def __init__(self, t, y, x, base):
        # Ensure that t, y, x are positive integers or throw an error
        if not(all(map(is_pos_int, (t, y, x)))):
            raise ValueError("The target dimensions (t, y, x) must all be positive integers.")

        self.ndim = 3
        self.size = t * y * x
        self.shape = (t, y, x)
        self.base = base

        match np.ndim(base):
            # 0D: Constant Input
            case 0:
                self.t_repetitions = t
                self.y_repetitions = y
                self.x_repetitions = x
                self.data = np.asarray([[[base]]])

            # 1D: Time-Series Input
            case 1:
                t0 = len(base)
                self.t_repetitions = divide_evenly(t, t0)
                self.y_repetitions = y
                self.x_repetitions = x
                # Expand (t0) -> (t0,1,1)
                self.data = np.expand_dims(base, axis=(1,2))

            # 2D: Spatial Input
            case 2:
                (y0, x0) = np.shape(base)
                self.t_repetitions = t
                self.y_repetitions = divide_evenly(y, y0)
                self.x_repetitions = divide_evenly(x, x0)
                # Expand (y0,x0) -> (1,y0,x0)
                self.data = np.expand_dims(base, axis=0)

            # 3D: Spatio-Temporal Input
            case 3:
                (t0, y0, x0) = np.shape(base)
                self.t_repetitions = divide_evenly(t, t0)
                self.y_repetitions = divide_evenly(y, y0)
                self.x_repetitions = divide_evenly(x, x0)
                self.data = np.asarray(base)

            # 4D+: Invalid Input
            case _:
                raise ValueError("Invalid input: base must have 0-3 dimensions.")


    def get(self, t, y, x):
        # Select value by spatio-temporal coordinate
        return self.data[t // self.t_repetitions,
                         y // self.y_repetitions,
                         x // self.x_repetitions]


    def getTimeSeries(self, y, x):
        # Select time series by spatial coordinate
        time = self.data[:,
                         y // self.y_repetitions,
                         x // self.x_repetitions]
        # Expand time dimension to (t)
        return maybe_repeat_array(time, (0, self.t_repetitions))


    def getSpatialPlane(self, t):
        # Select spatial plane by timestep
        space = self.data[t // self.t_repetitions]
        # Expand spatial dimensions to (y, x)
        return reduce(maybe_repeat_array,
                      ((0, self.y_repetitions),
                       (1, self.x_repetitions)),
                      space)


    def getSubCube(self, t_start, t_stop, y_start, y_stop, x_start, x_stop):
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
        return high_res_cube[t_start_idx:t_stop_idx,
                             y_start_idx:y_stop_idx,
                             x_start_idx:x_stop_idx]


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
        if hasattr(self, "cube"):
            return self.cube
        else:
            cube = self.__getFullyRealizedCube()
            if cache is True:
                self.cube = cube
            return cube


    def releaseFullyRealizedCube(self):
        delattr(self, "cube")
# space-time-cube-class ends here
