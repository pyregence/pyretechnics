cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy, vec_xyz, FuelModel, ProjectedVectors, FireBehaviorMin, FireBehaviorMax

cpdef float calc_flame_length(float fireline_intensity) noexcept

cpdef float calc_midflame_wind_speed(float wind_speed_20ft, float fuel_bed_depth, float canopy_height, float canopy_cover)

cpdef ProjectedVectors project_wind_and_slope_vectors_3d(
        float wind_speed, 
        float downwind_direction, 
        float slope,
        float upslope_direction
    ) noexcept

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

cpdef dict calc_surface_fire_behavior_in_direction(
    dict surface_fire_max,
    vec_xyz spread_direction,
    )

cpdef float calc_areal_heat_output(float spread_rate, float fireline_intensity) noexcept
