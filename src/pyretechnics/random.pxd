from pyretechnics.cy_types cimport pyidx

cdef class BufferedRandGen:
    cdef object numpy_rand
    
    cdef pyidx uniform_pos
    cdef float[:] uniform_buf

    cdef pyidx normal_pos
    cdef float[:] normal_buf

    cdef long[:] poisson16buf
    cdef pyidx poisson16pos
    cdef long[:] poisson1buf
    cdef pyidx poisson1pos
    cdef double[:] poisson_exp_buf
    cdef pyidx poisson_exp_pos

    cpdef float next_uniform(self)
    cpdef float next_normal(self)
    cpdef pyidx next_poisson(self, double M)

