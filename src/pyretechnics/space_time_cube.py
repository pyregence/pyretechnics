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
                # Repeat (t0) -> (t)              [FIXME: Inefficient memory usage]
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
                # Repeat (y0,x0) -> (y,x)         [FIXME: Inefficient memory usage]
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
                # Repeat (t0,y0,x0) -> (t,y,x)    [FIXME: Inefficient memory usage]
                repeated_array = reduce(SpaceTimeCube.maybe_repeat_array,
                                        ((0, t_repetitions),
                                         (1, y_repetitions),
                                         (2, x_repetitions)),
                                        base)
                self.data = repeated_array

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
# space-time-cube-class ends here
