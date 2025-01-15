cimport pyretechnics.cy_types
from pyretechnics.cy_types cimport vec_xyz, FireBehaviorMax, SpreadBehavior

#==================================================
# Cython functions to cimport into other modules
#==================================================

cdef float van_wagner_critical_fireline_intensity(
    float canopy_base_height,
    float foliar_moisture,
    ) noexcept

cdef float van_wagner_crowning_spread_rate_threshold(
    FireBehaviorMax surface_fire_max,
    float canopy_base_height,
    float foliar_moisture,
    )

cdef bint van_wagner_crown_fire_initiation(
    float surface_fireline_intensity,
    float canopy_cover,
    float canopy_base_height,
    float foliar_moisture,
    )

cdef FireBehaviorMax calc_crown_fire_behavior_max(
    float canopy_height,
    float canopy_base_height,
    float canopy_bulk_density,
    float heat_of_combustion,
    float estimated_fine_fuel_moisture,
    float wind_speed_10m,
    float upwind_direction,
    float slope,
    float aspect,
    float crown_max_lw_ratio=?,
    )

cdef SpreadBehavior calc_crown_fire_behavior_in_direction(
    FireBehaviorMax crown_fire_max,
    vec_xyz spread_direction,
    )

cdef SpreadBehavior calc_combined_fire_behavior(
    SpreadBehavior surface_fire_behavior,
    SpreadBehavior crown_fire_behavior,
    )
