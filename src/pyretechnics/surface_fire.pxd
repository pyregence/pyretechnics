cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xy, vec_xyz, FireBehaviorMax

cpdef float calc_flame_length(float fireline_intensity)

cpdef float calc_midflame_wind_speed(float wind_speed_20ft, float fuel_bed_depth, float canopy_height, float canopy_cover)

cpdef FireBehaviorMax calc_surface_fire_behavior_max(
    object surface_fire_min,
    float midflame_wind_speed,
    float upwind_direction,
    float slope,
    float aspect,
    bint use_wind_limit,
    object surface_lw_ratio_model
)
