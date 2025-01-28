from pyretechnics.cy_types cimport vec_xy, coord_yx, coord_tyx, SpreadBehavior
from pyretechnics.random cimport BufferedRandGen
from pyretechnics.space_time_cube cimport ISpaceTimeCube

cdef float expected_firebrand_production(
    SpreadBehavior fire_behavior,
    vec_xy elevation_gradient,
    float cell_horizontal_area,
    float firebrands_per_unit_heat=?,
    ) noexcept

cdef object spread_firebrands(
    ISpaceTimeCube fuel_model_cube,
    ISpaceTimeCube temperature_cube,
    ISpaceTimeCube fuel_moisture_dead_1hr_cube,
    unsigned char[:,:] fire_type_matrix,
    coord_yx sim_area_bounds,
    vec_xy spatial_resolution, 
    coord_tyx space_time_coordinate,
    float upwind_direction,
    float wind_speed_10m,
    float fireline_intensity,
    float flame_length,
    float time_of_arrival,
    BufferedRandGen random_generator,
    long num_firebrands,
    dict spot_config,
    )
