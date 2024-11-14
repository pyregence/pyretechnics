#==================================================
# C functions to cimport into other modules
#==================================================

cdef extern from "math.h":
     cdef double sqrt(double x)
     cdef double sin(double x)
     cdef double cos(double x)
