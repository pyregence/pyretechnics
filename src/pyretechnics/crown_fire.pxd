# [[file:../../org/pyretechnics.org::crown-fire-pxd][crown-fire-pxd]]
from pyretechnics.cy_types cimport vec_xyz, FireBehaviorMax, SpreadBehavior, CrownSpreadInfo

cdef float van_wagner_critical_fireline_intensity(
    float canopy_base_height,
    float foliar_moisture,
    ) noexcept

cdef float van_wagner_crowning_spread_rate(
    FireBehaviorMax surface_fire_max,
    float canopy_base_height,
    float foliar_moisture,
    ) noexcept

cpdef bint van_wagner_crown_fire_initiation(
    float surface_fireline_intensity,
    float canopy_cover,
    float canopy_base_height,
    float foliar_moisture,
    ) noexcept

cdef float cruz_active_crown_fire_spread_rate(
    float wind_speed_10m,
    float canopy_bulk_density,
    float estimated_fine_fuel_moisture,
    ) noexcept

cdef float van_wagner_critical_spread_rate(float canopy_bulk_density) noexcept

cdef float cruz_passive_crown_fire_spread_rate(
    float active_spread_rate,
    float critical_spread_rate,
    ) noexcept

cdef CrownSpreadInfo cruz_crown_fire_spread_info(
    float wind_speed_10m,
    float canopy_bulk_density,
    float estimated_fine_fuel_moisture,
    ) noexcept

cdef float calc_crown_fireline_intensity(
    float crown_spread_rate,
    float canopy_bulk_density,
    float canopy_height,
    float canopy_base_height,
    float heat_of_combustion,
    ) noexcept

cdef float crown_length_to_width_ratio(
    float wind_speed_10m,
    float max_length_to_width_ratio=?,
    ) noexcept

cdef float crown_fire_eccentricity(float length_to_width_ratio) noexcept

cpdef FireBehaviorMax calc_crown_fire_behavior_max(
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
    ) noexcept

cpdef SpreadBehavior calc_crown_fire_behavior_in_direction(
    FireBehaviorMax crown_fire_max,
    vec_xyz spread_direction,
    ) noexcept

cpdef SpreadBehavior calc_combined_fire_behavior(
    SpreadBehavior surface_fire_behavior,
    SpreadBehavior crown_fire_behavior,
    ) noexcept
# crown-fire-pxd ends here
