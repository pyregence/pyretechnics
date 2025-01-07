# [[file:../../org/pyretechnics.org::cy-types-pxd][cy-types-pxd]]
#==================================================
# Type aliases
#==================================================

ctypedef Py_ssize_t pyidx
ctypedef (float, float) vec_xy
ctypedef (float, float, float) vec_xyz
ctypedef (pyidx, pyidx) coord_yx
ctypedef (pyidx, pyidx, pyidx) coord_tyx

cdef struct FireBehaviorMax:
    int max_fire_type
    float max_spread_rate
    vec_xyz max_spread_direction
    float max_fireline_intensity
    float max_flame_length
    float length_to_width_ratio
    float eccentricity
    float critical_spread_rate
# cy-types-pxd ends here
