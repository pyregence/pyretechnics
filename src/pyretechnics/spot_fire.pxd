from pyretechnics.cy_types cimport vec_xy, SpreadBehavior
from pyretechnics.random cimport BufferedRandGen

cpdef float expct_firebrand_production(
    SpreadBehavior fire_behavior, 
    vec_xy elevation_gradient, 
    float cell_horizontal_area_m2, 
    float firebrands_per_unit_heat = ?
    )

cpdef object spread_firebrands(
    space_time_cubes, 
    output_matrices,
    cube_resolution, 
    space_time_coordinate,
    random_generator,
    expected_firebrand_count, 
    spot_config)