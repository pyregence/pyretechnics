from pyretechnics.cy_types cimport vec_xyz, fcatarr, fclaarr, FuelModel, ProjectedVectors, \
     FireBehaviorMin, FireBehaviorMax, SpreadBehavior

cdef float calc_surface_area_to_volume_ratio(fcatarr f_i, fclaarr f_ij, fclaarr sigma) noexcept

cdef float calc_packing_ratio(fclaarr w_o, fclaarr rho_p, float delta) noexcept

cdef float calc_optimum_packing_ratio(float sigma_prime) noexcept

cdef fcatarr calc_mineral_damping_coefficients(fclaarr f_ij, fclaarr S_e) noexcept

cdef fcatarr calc_moisture_damping_coefficients(fclaarr f_ij, fclaarr M_f, fclaarr M_x) noexcept

cdef fcatarr calc_low_heat_content(fclaarr f_ij, fclaarr h) noexcept

cdef fcatarr calc_net_fuel_loading(fclaarr g_ij, fclaarr w_o, fclaarr S_T) noexcept

cdef float calc_heat_per_unit_area(fcatarr eta_S_i, fcatarr eta_M_i, fcatarr h_i, fcatarr W_n_i) noexcept

cdef float calc_optimum_reaction_velocity(float sigma_prime, float beta, float beta_op) noexcept

cdef float calc_reaction_intensity(
    FuelModel moisturized_fuel_model,
    float sigma_prime,
    float beta,
    float beta_op,
    ) noexcept

cdef float calc_propagating_flux_ratio(float sigma_prime, float beta) noexcept

cdef float calc_heat_source(float I_R, float xi) noexcept

cdef float calc_ovendry_bulk_density(fclaarr w_o, float delta) noexcept

cdef fclaarr calc_effective_heating_number_distribution(fclaarr sigma) noexcept

cdef fclaarr calc_heat_of_preignition_distribution(fclaarr M_f) noexcept

cdef float calc_heat_sink(fcatarr f_i, fclaarr f_ij, float rho_b, fclaarr epsilon_ij, fclaarr Q_ig_ij) noexcept

cdef float calc_spread_rate(float heat_source, float heat_sink) noexcept

cdef float calc_residence_time(float sigma_prime) noexcept

cdef float calc_flame_depth(float spread_rate, float residence_time) noexcept

cdef float calc_fireline_intensity(float reaction_intensity, float flame_depth) noexcept

cdef float calc_flame_length(float fireline_intensity) noexcept

cdef float calc_areal_heat_output(float spread_rate, float fireline_intensity) noexcept

cdef float calc_max_effective_wind_speed(float reaction_intensity) noexcept

cdef float get_phi_S(FireBehaviorMin sfmin, float slope) noexcept

cdef float get_phi_W(FireBehaviorMin sfmin, float midflame_wind_speed) noexcept

cdef float get_wind_speed(FireBehaviorMin sfmin, float phi_W) noexcept

cdef FireBehaviorMin make_surface_fire_min(
    float base_spread_rate,
    float base_fireline_intensity,
    float max_effective_wind_speed,
    float B,
    float C,
    float F,
    float beta,
    ) noexcept

cpdef FireBehaviorMin calc_surface_fire_behavior_no_wind_no_slope(
    FuelModel moisturized_fuel_model,
    float spread_rate_adjustment=?,
    ) noexcept

cdef float calc_wind_adjustment_factor(float fuel_bed_depth, float canopy_height, float canopy_cover) noexcept

cpdef float calc_midflame_wind_speed(
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

cdef vec_xyz get_phi_E(vec_xyz wind_vector_3d, vec_xyz slope_vector_3d, float phi_W, float phi_S) noexcept

cdef float surface_length_to_width_ratio(float effective_wind_speed, str model=?)

cdef float surface_fire_eccentricity(float length_to_width_ratio) noexcept

cdef (float, float) maybe_limit_wind_speed(
    bint use_wind_limit,
    float max_wind_speed,
    FireBehaviorMin sfmin,
    float phi_E_magnitude,
    ) noexcept

cpdef FireBehaviorMax calc_surface_fire_behavior_max(
    FireBehaviorMin sfmin,
    float midflame_wind_speed,
    float upwind_direction,
    float slope,
    float aspect,
    bint use_wind_limit=?,
    str surface_lw_ratio_model=?,
    ) noexcept

cpdef SpreadBehavior calc_surface_fire_behavior_in_direction(
    FireBehaviorMax surface_fire_max,
    vec_xyz spread_direction,
    ) noexcept
