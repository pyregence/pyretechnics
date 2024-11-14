cimport pyretechnics.types
from pyretechnics.types cimport vec_xy

#==================================================
# Cython functions to cimport into other modules
#==================================================

cpdef float dot_2d(vec_xy vector1, vec_xy vector2)
cpdef float vector_magnitude_2d(vec_xy vector)
