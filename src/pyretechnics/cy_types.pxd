# [[file:../../org/pyretechnics.org::cy-types-pxd][cy-types-pxd]]
#==================================================
# Compile-time-defined type aliases
#==================================================

ctypedef Py_ssize_t pyidx
ctypedef (float, float) vec_xy
ctypedef (float, float, float) vec_xyz
ctypedef (pyidx, pyidx) coord_yx
ctypedef (pyidx, pyidx, pyidx) coord_tyx
ctypedef (float, float) fcatarr
ctypedef (float, float, float, float, float, float) fclaarr

cdef struct FuelModel:
    int number
    float delta
    fclaarr M_x
    fclaarr M_f
    fclaarr w_o
    fclaarr sigma
    fclaarr h
    fclaarr rho_p
    fclaarr S_T
    fclaarr S_e
    bint dynamic
    bint burnable
    fclaarr exp_A_sigma
    fclaarr firemod_size_classes
    fclaarr f_ij
    fcatarr f_i
    fclaarr g_ij

cdef struct FireBehaviorMax:
    int max_fire_type
    float max_spread_rate
    vec_xyz max_spread_direction
    float max_fireline_intensity
    float max_flame_length
    float length_to_width_ratio
    float eccentricity
    float critical_spread_rate

cdef struct SpreadBehavior:
    float dphi_dt
    int fire_type
    float spread_rate
    vec_xyz spread_direction
    float fireline_intensity
    float flame_length
# cy-types-pxd ends here
