# [[file:../../org/pyretechnics.org::random-pxd][random-pxd]]
from pyretechnics.cy_types cimport pyidx

cdef class BufferedRandGen:
    cdef object numpy_rand     # TODO: object -> Generator
    cdef float[::1] uniform_buf
    cdef pyidx uniform_pos
    cdef float[::1] normal_buf
    cdef pyidx normal_pos
    cdef long[::1] poisson16_buf
    cdef pyidx poisson16_pos
    cdef long[::1] poisson1_buf
    cdef pyidx poisson1_pos
    cdef double[::1] poisson_exp_buf
    cdef pyidx poisson_exp_pos
    cdef long next_poisson(BufferedRandGen self, double M) noexcept
    cdef float next_uniform(BufferedRandGen self) noexcept
    cdef float next_normal(BufferedRandGen self) noexcept
# random-pxd ends here
