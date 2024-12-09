from pyretechnics.cy_types cimport pyidx

cdef class BufferedRandGen:
    cdef object numpy_rand
    
    cdef pyidx uniform_pos
    cdef float[:] uniform_buf

    cdef pyidx normal_pos
    cdef float[:] normal_buf

    cpdef float next_uniform(self)
    cpdef float next_normal(self)

