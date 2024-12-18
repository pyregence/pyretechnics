from pyretechnics.cy_types cimport pyidx

cdef class ISpaceTimeCube:
    cpdef float get(self, pyidx t, pyidx y, pyidx x)
    cdef float get_float(self, pyidx t, pyidx y, pyidx x) # FIXME return double


cdef class SpaceTimeCube(ISpaceTimeCube):
    cdef public int ndim
    cdef public unsigned long size
    cdef public tuple shape
    cdef public object base
    cdef public int t_repetitions
    cdef public int y_repetitions
    cdef public int x_repetitions
    cdef public float[:,:,:] data

