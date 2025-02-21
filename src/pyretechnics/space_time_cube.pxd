from pyretechnics.cy_types cimport pyidx

cdef bint is_pos_int(object x) noexcept
cdef int divide_evenly(int dividend, int divisor)

# def to_positive_index_range(index_range, axis_length)
# def maybe_repeat_array(array, axis_repetitions)

cdef class ISpaceTimeCube:
    cpdef float get(ISpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept

cdef class SpaceTimeCube(ISpaceTimeCube):
    cdef public int ndim
    cdef public unsigned long size
    cdef public (int, int, int) shape
    cdef public object base
    cdef public int t_repetitions
    cdef public int y_repetitions
    cdef public int x_repetitions
    cdef public float[:,:,:] data
    cpdef float get(SpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept
    # def getTimeSeries(self, t_range, y, x)
    # def getSpatialPlane(self, t, y_range, x_range)
    # def getSubcube(self, t_range, y_range, x_range)
    # def getFullyRealizedCube(self, cache=False)

cdef class LazySpaceTimeCube(ISpaceTimeCube):
    cdef public int ndim
    cdef public unsigned long size
    cdef public (int, int, int) shape
    cdef public (int, int, int) subcube_shape
    cdef public (int, int, int) cache_shape
    cdef public object[:,:,:] cache
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
