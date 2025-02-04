from pyretechnics.cy_types cimport pyidx

cdef class ISpaceTimeCube:
    cdef float get(ISpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept

cdef class SpaceTimeCube(ISpaceTimeCube):
    cdef public int ndim
    cdef public unsigned long size
    cdef public (int, int, int) shape
    cdef public object base
    cdef public int t_repetitions
    cdef public int y_repetitions
    cdef public int x_repetitions
    cdef public float[:,:,:] data
    cdef float get(SpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept

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
    cdef float get(LazySpaceTimeCube self, pyidx t, pyidx y, pyidx x) noexcept
