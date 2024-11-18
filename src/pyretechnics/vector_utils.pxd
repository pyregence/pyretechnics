# [[file:../../org/pyretechnics.org::vector-utils-pxd][vector-utils-pxd]]
cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy, vec_xyz

#==================================================
# Cython functions to cimport into other modules
#==================================================

cpdef double vector_magnitude(double[:] vector)
cpdef float dot_2d(vec_xy vector1, vec_xy vector2)
cpdef float dot_3d(vec_xyz vector1, vec_xyz vector2)
cpdef float vector_magnitude_2d(vec_xy vector)
cpdef float vector_magnitude_3d(vec_xyz vector)
cpdef vec_xy as_unit_vector_2d(vec_xy vector)
cpdef vec_xyz as_unit_vector_3d(vec_xyz vector)
cpdef vec_xyz to_slope_plane(vec_xy vector_2d, vec_xy elevation_gradient)
cpdef double[:] to_horizontal_plane(double[:] vector_3d)
cpdef float spread_direction_vector_to_angle(double[:] vector_3d)
cpdef vec_xyz get_slope_normal_vector(vec_xy elevation_gradient)
cpdef double[:] rotate_on_sloped_plane(double[:] vector, float theta, float slope, float aspect)
# vector-utils-pxd ends here
