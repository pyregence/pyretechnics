cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy, vec_xyz, FuelModel, FireBehaviorMax

cpdef float calc_flame_length(float fireline_intensity)

cpdef float calc_midflame_wind_speed(float wind_speed_20ft, float fuel_bed_depth, float canopy_height, float canopy_cover)


cdef struct ProjectedVectors:
    vec_xyz wind_vector_3d
    vec_xyz slope_vector_3d


cpdef ProjectedVectors project_wind_and_slope_vectors_3d(
        float wind_speed, 
        float downwind_direction, 
        float slope,
        float upslope_direction
    )


cdef struct FireBehaviorMin:
    float base_spread_rate
    float base_fireline_intensity
    float max_effective_wind_speed
    float _phiS_G
    float _phiW_scalr
    float _phiW_expnt
    float _ws_scalr
    float _ws_expnt


cdef FireBehaviorMin calc_surface_fire_behavior_no_wind_no_slope(
    FuelModel moisturized_fuel_model, 
    float spread_rate_adjustment)


cpdef FireBehaviorMax calc_surface_fire_behavior_max(
    FireBehaviorMin surface_fire_min,
    float midflame_wind_speed,
    float upwind_direction,
    float slope,
    float aspect,
    bint use_wind_limit,
    object surface_lw_ratio_model
    )


