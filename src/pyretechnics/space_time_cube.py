# [[file:../../org/pyretechnics.org::space-time-cube-class][space-time-cube-class]]
from functools import reduce
import numpy as np

class SpaceTimeCube:
    """
    TODO: Add docstring.
    """
    def __init__(self, base, t, y, x):
        # Ensure that t, y, x are positive integers or throw an error
        if not(all(map(SpaceTimeCube.is_pos_int, (t, y, x)))):
            raise ValueError("The target dimensions (t, y, x) must all be positive integers.")

        self.ndim = 3
        self.size = t * y * x
        self.shape = (t, y, x)
        self.base = base

        match np.ndim(base):
            # 0D: Constant Input
            case 0:
                self.data = np.broadcast_to([[[base]]], self.shape)

            # 1D: Time-Series Input
            case 1:
                t0 = len(base)
                t_repetitions = SpaceTimeCube.divide_evenly(t, t0)
                # Repeat (t0) -> (t)              [NOTE: Inefficient memory usage]
                repeated_array = SpaceTimeCube.maybe_repeat_array(base, (0, t_repetitions))
                # Expand (t) -> (t,1,1)
                expanded_array = np.expand_dims(repeated_array, axis=(1,2))
                # Broadcast (t,1,1) -> (t,y,x)
                self.data = np.broadcast_to(expanded_array, self.shape)

            # 2D: Spatial Input
            case 2:
                (y0, x0) = np.shape(base)
                y_repetitions = SpaceTimeCube.divide_evenly(y, y0)
                x_repetitions = SpaceTimeCube.divide_evenly(x, x0)
                # Repeat (y0,x0) -> (y,x)         [NOTE: Inefficient memory usage]
                repeated_array = reduce(SpaceTimeCube.maybe_repeat_array,
                                        ((0, y_repetitions),
                                         (1, x_repetitions)),
                                        base)
                # Expand (y,x) -> (1,y,x)
                expanded_array = np.expand_dims(repeated_array, axis=0)
                # Broadcast (1,y,x) -> (t,y,x)
                self.data = np.broadcast_to(expanded_array, self.shape)

            # 3D: Spatio-Temporal Input
            case 3:
                (t0, y0, x0) = np.shape(base)
                t_repetitions = SpaceTimeCube.divide_evenly(t, t0)
                y_repetitions = SpaceTimeCube.divide_evenly(y, y0)
                x_repetitions = SpaceTimeCube.divide_evenly(x, x0)
                # Repeat (t0,y0,x0) -> (t,y,x)    [NOTE: Inefficient memory usage]
                repeated_array = reduce(SpaceTimeCube.maybe_repeat_array,
                                        ((0, t_repetitions),
                                         (1, y_repetitions),
                                         (2, x_repetitions)),
                                        base)
                self.data = np.asarray(repeated_array)

            # 4D+: Invalid Input
            case _:
                raise ValueError("Invalid input: base must have 0-3 dimensions.")

    @staticmethod
    def is_pos_int(x):
        return isinstance(x, int) and x > 0

    @staticmethod
    def divide_evenly(dividend, divisor):
        (quotient, remainder) = divmod(dividend, divisor)
        if remainder == 0:
            return quotient
        else:
            raise ValueError(str(dividend) + " must be an exact multiple of " + str(divisor) + ".")

    @staticmethod
    def maybe_repeat_array(array, axis_repetitions):
        """
        TODO: Add docstring.
        """
        (axis, repetitions) = axis_repetitions
        return np.repeat(array, repetitions, axis) if (repetitions > 1) else array

    def getSpatialPlane(self, t):
        """
        TODO: Add docstring.
        """
        return self.data[t]

    def getTimeSeries(self, y, x):
        """
        TODO: Add docstring.
        """
        return self.data[:,y,x]

    def get(self, t, y, x):
        """
        TODO: Add docstring.
        """
        return self.data[t,y,x]

    def getSubCube(self, t_start, t_stop, y_start, y_stop, x_start, x_stop):
        """
        TODO: Add docstring.
        """
        return self.data[t_start:t_stop, y_start:y_stop, x_start:x_stop]
# space-time-cube-class ends here
# [[file:../../org/pyretechnics.org::space-time-cube2-class][space-time-cube2-class]]
from functools import reduce
import numpy as np

class SpaceTimeCube2:
    """
    TODO: Add docstring.
    """
    def __init__(self, base, t, y, x):
        # Ensure that t, y, x are positive integers or throw an error
        if not(all(map(SpaceTimeCube2.is_pos_int, (t, y, x)))):
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
                self.t_repetitions = SpaceTimeCube2.divide_evenly(t, t0)
                self.y_repetitions = y
                self.x_repetitions = x
                self.data = np.expand_dims(base, axis=(1,2))

            # 2D: Spatial Input
            case 2:
                (y0, x0) = np.shape(base)
                self.t_repetitions = t
                self.y_repetitions = SpaceTimeCube2.divide_evenly(y, y0)
                self.x_repetitions = SpaceTimeCube2.divide_evenly(x, x0)
                self.data = np.expand_dims(base, axis=0)

            # 3D: Spatio-Temporal Input
            case 3:
                (t0, y0, x0) = np.shape(base)
                self.t_repetitions = SpaceTimeCube2.divide_evenly(t, t0)
                self.y_repetitions = SpaceTimeCube2.divide_evenly(y, y0)
                self.x_repetitions = SpaceTimeCube2.divide_evenly(x, x0)
                self.data = np.asarray(base)

            # 4D+: Invalid Input
            case _:
                raise ValueError("Invalid input: base must have 0-3 dimensions.")

    @staticmethod
    def is_pos_int(x):
        return isinstance(x, int) and x > 0

    @staticmethod
    def divide_evenly(dividend, divisor):
        (quotient, remainder) = divmod(dividend, divisor)
        if remainder == 0:
            return quotient
        else:
            raise ValueError(str(dividend) + " must be an exact multiple of " + str(divisor) + ".")

    @staticmethod
    def maybe_repeat_array(array, axis_repetitions):
        """
        TODO: Add docstring.
        """
        (axis, repetitions) = axis_repetitions
        return np.repeat(array, repetitions, axis) if (repetitions > 1) else array

    def getSpatialPlane(self, t):
        """
        TODO: Add docstring.
        """
        space = self.data[t // self.t_repetitions]
        # [NOTE: CPU and Memory Inefficient]
        return reduce(SpaceTimeCube2.maybe_repeat_array,
                      ((0, self.y_repetitions),
                       (1, self.x_repetitions)),
                      space)

    def getTimeSeries(self, y, x):
        """
        TODO: Add docstring.
        """
        time = self.data[:,
                         y // self.y_repetitions,
                         x // self.x_repetitions]
        # [NOTE: CPU and Memory Inefficient]
        return SpaceTimeCube2.maybe_repeat_array(time, (0, self.t_repetitions))

    def get(self, t, y, x):
        """
        TODO: Add docstring.
        """
        return self.data[t // self.t_repetitions,
                         y // self.y_repetitions,
                         x // self.x_repetitions]

    def getSubCube(self, t_start, t_stop, y_start, y_stop, x_start, x_stop):
        """
        TODO: Add docstring.
        """
        t_start_chunk = t_start // self.t_repetitions
        t_stop_chunk  = t_stop  // self.t_repetitions
        y_start_chunk = y_start // self.y_repetitions
        y_stop_chunk  = y_stop  // self.y_repetitions
        x_start_chunk = x_start // self.x_repetitions
        x_stop_chunk  = x_stop  // self.x_repetitions
        low_res_cube  = self.data[t_start_chunk:(t_stop_chunk + 1),
                                  y_start_chunk:(y_stop_chunk + 1),
                                  x_start_chunk:(x_stop_chunk + 1)]
        # [NOTE: CPU and Memory Inefficient]
        high_res_cube = reduce(SpaceTimeCube2.maybe_repeat_array,
                               ((0, self.t_repetitions),
                                (1, self.y_repetitions),
                                (2, self.x_repetitions)),
                               low_res_cube)
        t_chunk_origin = t_start_chunk * self.t_repetitions
        y_chunk_origin = y_start_chunk * self.y_repetitions
        x_chunk_origin = x_start_chunk * self.x_repetitions
        t_start_idx    = t_start - t_chunk_origin
        t_stop_idx     = t_stop  - t_chunk_origin
        y_start_idx    = y_start - y_chunk_origin
        y_stop_idx     = y_stop  - y_chunk_origin
        x_start_idx    = x_start - x_chunk_origin
        x_stop_idx     = x_stop  - x_chunk_origin
        return high_res_cube[t_start_idx:t_stop_idx,
                             y_start_idx:y_stop_idx,
                             x_start_idx:x_stop_idx]
# space-time-cube2-class ends here
