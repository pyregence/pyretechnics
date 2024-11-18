# [[file:../../org/pyretechnics.org::conversion-pxd][conversion-pxd]]
cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy

#==================================================
# Cython functions to cimport into other modules
#==================================================

cpdef float rad_to_deg(float radians)
cpdef float deg_to_rad(float degrees)
cpdef vec_xy cartesian_to_azimuthal(float x, float y)
cpdef vec_xy azimuthal_to_cartesian(float r, float azimuth)
cpdef float opposite_direction(float theta)
# conversion-pxd ends here
