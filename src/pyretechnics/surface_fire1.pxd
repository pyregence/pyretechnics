cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy, vec_xyz, FuelModel, FireBehaviorMax

cpdef float calc_flame_length(float fireline_intensity)

cpdef float calc_midflame_wind_speed(float wind_speed_20ft, float fuel_bed_depth, float canopy_height, float canopy_cover)

cdef class FireBehaviorMin:
    cdef public float base_spread_rate
    cdef public float base_fireline_intensity
    cdef public float max_effective_wind_speed
    cdef float _phiS_G
    cdef float _phiW_scalr
    cdef float _phiW_expnt
    cdef float _ws_scalr
    cdef float _ws_expnt
    cpdef float get_phi_W(self, float midflame_wind_speed)
    cpdef float get_phi_S(self, float slope)
    cpdef float get_wind_speed(self, float phi_W)


cdef FireBehaviorMin calc_surface_fire_behavior_no_wind_no_slope(
    #FuelModel moisturized_fuel_model, 
    object moisturized_fuel_model, 
    float spread_rate_adjustment)


cpdef FireBehaviorMax calc_surface_fire_behavior_max(
    object surface_fire_min,
    float midflame_wind_speed,
    float upwind_direction,
    float slope,
    float aspect,
    bint use_wind_limit,
    object surface_lw_ratio_model
    )


