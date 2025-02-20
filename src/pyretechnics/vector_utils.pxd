# [[file:../../org/pyretechnics.org::vector-utils-pxd][vector-utils-pxd]]
from pyretechnics.cy_types cimport vec_xy, vec_xyz

cdef float dot_2d(vec_xy vector1, vec_xy vector2) noexcept
cdef float dot_3d(vec_xyz vector1, vec_xyz vector2) noexcept
cdef vec_xy scale_2d(float scalar, vec_xy vector) noexcept
cdef vec_xyz scale_3d(float scalar, vec_xyz vector) noexcept
cdef vec_xy add_2d(vec_xy vector1, vec_xy vector2) noexcept
cdef vec_xyz add_3d(vec_xyz vector1, vec_xyz vector2) noexcept
cdef float vector_magnitude_2d(vec_xy vector) noexcept
cdef float vector_magnitude_3d(vec_xyz vector) noexcept
cdef vec_xy as_unit_vector_2d(vec_xy vector) noexcept
cdef vec_xyz as_unit_vector_3d(vec_xyz vector) noexcept
cdef vec_xyz to_slope_plane(vec_xy vector_2d, vec_xy elevation_gradient) noexcept
cdef vec_xy to_horizontal_plane(vec_xyz vector_3d) noexcept
cpdef float spread_direction_vector_to_angle(vec_xyz vector_3d) noexcept
cdef vec_xyz get_slope_normal_vector(vec_xy elevation_gradient) noexcept
cdef vec_xyz cross_3d(vec_xyz vector1, vec_xyz vector2) noexcept
cdef vec_xyz rotate_on_sloped_plane(vec_xyz vector, float theta, float slope, float aspect) noexcept
# vector-utils-pxd ends here
