from pyretechnics.cy_types cimport coord_tyx, SpreadBehavior

cpdef SpreadBehavior burn_cell_as_head_fire(
    dict space_time_cubes,
    coord_tyx space_time_coordinate,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    float crown_max_lw_ratio=?,
    ) noexcept

cpdef SpreadBehavior burn_cell_toward_azimuth(
    dict space_time_cubes,
    coord_tyx space_time_coordinate,
    float azimuth,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    float crown_max_lw_ratio=?,
    ) noexcept
