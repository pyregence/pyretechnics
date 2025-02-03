from pyretechnics.cy_types cimport vec_xyz, FuelModel, ProjectedVectors, \
     FireBehaviorMin, FireBehaviorMax, SpreadBehavior

cdef float calc_flame_length(float fireline_intensity) noexcept

cdef float calc_midflame_wind_speed(
    float wind_speed_20ft,
    float fuel_bed_depth,
    float canopy_height,
    float canopy_cover,
    ) noexcept

cdef ProjectedVectors project_wind_and_slope_vectors_3d(
    float wind_speed,
    float downwind_direction,
    float slope,
    float upslope_direction,
    ) noexcept

cdef FireBehaviorMin calc_surface_fire_behavior_no_wind_no_slope(
    FuelModel moisturized_fuel_model,
    float spread_rate_adjustment,
    ) noexcept

cdef FireBehaviorMax calc_surface_fire_behavior_max(
    FireBehaviorMin surface_fire_min,
    float midflame_wind_speed,
    float upwind_direction,
    float slope,
    float aspect,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    ) noexcept

cdef SpreadBehavior calc_surface_fire_behavior_in_direction(
    FireBehaviorMax surface_fire_max,
    vec_xyz spread_direction,
    ) noexcept

cdef float calc_areal_heat_output(float spread_rate, float fireline_intensity) noexcept
