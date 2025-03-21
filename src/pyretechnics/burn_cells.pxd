# [[file:../../org/pyretechnics.org::burn-cells-pxd][burn-cells-pxd]]
from pyretechnics.cy_types cimport pyidx, coord_tyx

cpdef dict burn_cell_as_head_fire(
    dict space_time_cubes,
    coord_tyx space_time_coordinate,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    float crown_max_lw_ratio=?,
    )

cpdef dict burn_all_cells_as_head_fire(
    dict space_time_cubes,
    pyidx t,
    tuple y_range=?,
    tuple x_range=?,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    float crown_max_lw_ratio=?,
    )

cpdef dict burn_cell_toward_azimuth(
    dict space_time_cubes,
    coord_tyx space_time_coordinate,
    float azimuth,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    float crown_max_lw_ratio=?,
    )

cpdef dict burn_all_cells_toward_azimuth(
    dict space_time_cubes,
    float azimuth,
    pyidx t,
    tuple y_range=?,
    tuple x_range=?,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    float crown_max_lw_ratio=?,
    )
# burn-cells-pxd ends here
