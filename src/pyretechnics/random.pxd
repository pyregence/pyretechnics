from pyretechnics.cy_types cimport pyidx

cdef class BufferedRandGen:
    cdef object numpy_rand     # TODO: object -> Generator
    cdef float[:] uniform_buf
    cdef pyidx uniform_pos
    cdef float[:] normal_buf
    cdef pyidx normal_pos
    cdef long[:] poisson16_buf
    cdef pyidx poisson16_pos
    cdef long[:] poisson1_buf
    cdef pyidx poisson1_pos
    cdef double[:] poisson_exp_buf
    cdef pyidx poisson_exp_pos
    cdef pyidx next_poisson(BufferedRandGen self, double M) noexcept
    cdef float next_uniform(BufferedRandGen self) noexcept
    cdef float next_normal(BufferedRandGen self) noexcept
