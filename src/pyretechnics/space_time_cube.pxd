# [[file:../../org/pyretechnics.org::space-time-cube-pxd][space-time-cube-pxd]]
cimport numpy as np
from pyretechnics.cy_types cimport pyidx

cdef bint is_pos_int(object x) noexcept
cdef int divide_evenly(int dividend, int divisor)
cpdef (pyidx, pyidx) to_positive_index_range(object index_range, pyidx axis_length) noexcept
cdef np.ndarray maybe_repeat_array(np.ndarray array, (pyidx, int) axis_repetitions)

cdef class ISpaceTimeCube:
    cdef float get(ISpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept

cdef class SpaceTimeCube(ISpaceTimeCube):
    cdef public int ndim
    cdef public unsigned long long size
    cdef public (int, int, int) shape
    cdef public object base
    cdef public int t_repetitions
    cdef public int y_repetitions
    cdef public int x_repetitions
    # NOTE: We use const (read-only MemoryView) so that Cython will accept read-only arrays, which is required for shared-memory parallelism.
    cdef public const float[:,:,::1] data
    cpdef float get(SpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept
    cpdef np.ndarray getTimeSeries(SpaceTimeCube self, object t_range, pyidx y, pyidx x)
    cpdef np.ndarray getSpatialPlane(SpaceTimeCube self, pyidx t, object y_range, object x_range)
    cpdef np.ndarray getSubcube(SpaceTimeCube self, object t_range, object y_range, object x_range)
    # def getFullyRealizedCube(self, cache=False)
    # def releaseFullyRealizedCube(self)

cdef class LazySpaceTimeCube(ISpaceTimeCube):
    cdef public int ndim
    cdef public unsigned long long size
    cdef public (int, int, int) shape
    cdef public (int, int, int) subcube_shape
    cdef public (int, int, int) cache_shape
    cdef public np.ndarray cache
    cdef public object load_subcube
    cdef SpaceTimeCube __getOrLoadSubcube(
        LazySpaceTimeCube self,
        pyidx cache_t,
        pyidx cache_y,
        pyidx cache_x,
        )
    cpdef float get(LazySpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept
    # def getTimeSeries(self, t_range, y, x)
    # def getSpatialPlane(self, t, y_range, x_range)
    # def getSubcube(self, t_range, y_range, x_range)
    # def getFullyRealizedCube(self, cache=False)
    # def releaseFullyRealizedCube(self)
# space-time-cube-pxd ends here
