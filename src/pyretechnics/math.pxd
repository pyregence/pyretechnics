# [[file:../../org/pyretechnics.org::math-pxd][math-pxd]]
#==================================================
# C functions to cimport into other modules
#==================================================

cdef extern from "math.h":
     cdef double sqrt(double x)
     cdef double sin(double x)
     cdef double cos(double x)
     cdef double log(double x)
     cdef double exp(double x)
     cdef double pow(double x, double x)
     cdef double tan(double x)
     cdef double atan(double x)
     cdef double atan2(double y, double x)
# math-pxd ends here
