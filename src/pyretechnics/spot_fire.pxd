# [[file:../../org/pyretechnics.org::spot-fire-pxd][spot-fire-pxd]]
from pyretechnics.cy_types cimport pyidx, vec_xy, coord_yx, coord_tyx, SpreadBehavior, SpotConfig, JumpDistribution
from pyretechnics.random cimport BufferedRandGen
from pyretechnics.space_time_cube cimport ISpaceTimeCube

cdef float expected_firebrand_production(
    SpreadBehavior fire_behavior,
    vec_xy elevation_gradient,
    float cell_horizontal_area,
    float firebrands_per_unit_heat=?,
    ) noexcept

cdef float delta_to_grid_dx(float cos_wdir, float sin_wdir, float delta_x, float delta_y) noexcept
cdef float delta_to_grid_dy(float cos_wdir, float sin_wdir, float delta_x, float delta_y) noexcept

cdef int distance_to_n_cells(float distance, float cell_size) noexcept

cdef float resolve_exp_delta_x(SpotConfig spot_config, float fireline_intensity, float wind_speed_20ft) noexcept
cdef float resolve_var_delta_x(SpotConfig spot_config, float exp_delta_x) noexcept

cdef float lognormal_mu_from_moments(float mean, float variance) noexcept
cdef float lognormal_sigma_from_moments(float mean, float variance) noexcept

cdef (float, float) resolve_lognormal_params(
    SpotConfig spot_config,
    float fireline_intensity,
    float wind_speed_20ft
    ) noexcept

cdef double sigma_y_scalar_m = 0.92 * 0.47 / (0.88 * 0.88)

cdef float himoto_resolve_default_sigma_y_from_lognormal_params(float mu_x, float sigma_x) noexcept
cdef float himoto_resolve_default_sigma_y(
    SpotConfig spot_config,
    float fireline_intensity,
    float wind_speed_20ft
    ) noexcept

cdef float resolve_crosswind_distance_stdev(
    SpotConfig spot_config,
    float fireline_intensity,
    float wind_speed_20ft
    ) noexcept

cdef float sample_normal(BufferedRandGen rng, float mu, float sd) noexcept
cdef float sample_lognormal(BufferedRandGen rng, float mu, float sd) noexcept

cdef JumpDistribution resolve_JumpDistribution(SpotConfig spot_config, float fireline_intensity, float wind_speed_20ft)

cdef float sample_downwind_jump(JumpDistribution jd, BufferedRandGen random_generator) noexcept
cdef float sample_crosswind_jump(JumpDistribution jd, BufferedRandGen random_generator) noexcept

cdef float heat_of_preignition(float temperature, float fine_fuel_moisture) noexcept
cpdef float schroeder_ignition_probability(float temperature, float fine_fuel_moisture) noexcept
cpdef float firebrand_flight_survival_probability(float spotting_distance, float decay_distance) noexcept

cdef float albini_firebrand_maximum_height(float firebrand_diameter) noexcept
cdef float albini_t_max(float flame_length) noexcept
cdef float spot_ignition_time(float time_of_arrival, float flame_length) noexcept

cdef bint is_in_bounds(pyidx y, pyidx x, pyidx rows, pyidx cols) noexcept
cdef bint is_burnable_cell(ISpaceTimeCube fuel_model_cube, pyidx t, pyidx y, pyidx x) noexcept

cdef coord_yx cast_firebrand(
    BufferedRandGen rng,
    ISpaceTimeCube fuel_model_cube,
    ISpaceTimeCube temperature_cube,
    ISpaceTimeCube fuel_moisture_dead_1hr_cube,
    unsigned char[:,::1] fire_type_matrix,
    pyidx rows,
    pyidx cols,
    float cell_height,
    float cell_width,
    pyidx source_t,
    pyidx source_y,
    pyidx source_x,
    float decay_distance,
    float cos_wdir,
    float sin_wdir,
    JumpDistribution jd,
    ) noexcept

cdef tuple spread_firebrands(
    ISpaceTimeCube fuel_model_cube,
    ISpaceTimeCube temperature_cube,
    ISpaceTimeCube fuel_moisture_dead_1hr_cube,
    unsigned char[:,::1] fire_type_matrix,
    coord_yx sim_area_bounds,
    float cell_height,
    float cell_width,
    coord_tyx space_time_coordinate,
    float wind_speed_10m,
    float upwind_direction,
    float fireline_intensity,
    float flame_length,
    float time_of_arrival,
    BufferedRandGen random_generator,
    long long num_firebrands,
    SpotConfig spot_config,
    )
# spot-fire-pxd ends here
