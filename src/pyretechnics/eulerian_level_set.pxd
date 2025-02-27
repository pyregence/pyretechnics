# from pyretechnics.cy_types cimport pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, \
#      FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior, SpotConfig, \
#      PartialedEllWavelet, CellInputs, EllipticalInfo, Pass1CellOutput
# from pyretechnics.space_time_cube cimport ISpaceTimeCube
# from pyretechnics.narrow_band_tracking cimport NarrowBandTracker
# from pyretechnics.random cimport BufferedRandGen
# import numpy as np
#
# cdef float calc_dphi_dx_approx(float[:,:] phi_matrix, float dx, pyidx x, pyidx y) noexcept
# cdef float calc_dphi_dy_approx(float[:,:] phi_matrix, float dy, pyidx x, pyidx y) noexcept
# cdef vec_xy calc_phi_gradient_approx(float[:,:] phi_matrix, float dx, float dy, pyidx x, pyidx y) noexcept
# cdef vec_xy calc_phi_normal_vector(float[:,:] phi_matrix, float dx, float dy, pyidx x, pyidx y) noexcept
# cdef float calc_phi_normal_azimuth(vec_xy phi_normal_vector) noexcept
# cdef float half_superbee_dphi_up(float dphi_up, float dphi_loc) noexcept
# cdef float calc_dphi_flim_x(float p00, float pw2, float pw1, float pe1, float pe2) noexcept
# cdef float calc_dphi_flim_y(float p00, float ps2, float ps1, float pn1, float pn2) noexcept
# cdef int fire_type_unburned      = 0
# cdef int fire_type_surface       = 1
# cdef int fire_type_crown_passive = 2
# cdef int fire_type_crown_active  = 3
# cdef vec_xy calc_elevation_gradient(float slope, float aspect) noexcept
# cdef vec_xyz calc_phi_gradient_on_slope(vec_xy phi_gradient_xy, vec_xy elevation_gradient) noexcept
# cdef SpreadBehavior calc_fireline_normal_behavior(FireBehaviorMax fire_behavior_max, vec_xyz phi_gradient) noexcept
# cdef class SpreadInputs:
#     cdef pyidx rows
#     cdef pyidx cols
#     cdef float band_duration
#     cdef vec_xy spatial_resolution
#     cdef ISpaceTimeCube slope
#     cdef ISpaceTimeCube aspect
#     cdef ISpaceTimeCube fuel_model
#     cdef ISpaceTimeCube canopy_cover
#     cdef ISpaceTimeCube canopy_height
#     cdef ISpaceTimeCube canopy_base_height
#     cdef ISpaceTimeCube canopy_bulk_density
#     cdef ISpaceTimeCube wind_speed_10m
#     cdef ISpaceTimeCube upwind_direction
#     cdef ISpaceTimeCube fuel_moisture_dead_1hr
#     cdef ISpaceTimeCube fuel_moisture_dead_10hr
#     cdef ISpaceTimeCube fuel_moisture_dead_100hr
#     cdef ISpaceTimeCube fuel_moisture_live_herbaceous
#     cdef ISpaceTimeCube fuel_moisture_live_woody
#     cdef ISpaceTimeCube foliar_moisture
#     cdef ISpaceTimeCube temperature
#     cdef ISpaceTimeCube fuel_spread_adjustment
#     cdef ISpaceTimeCube weather_spread_adjustment
#     cdef FuelModel* fuel_models_arr
#     cdef void __init_fuel_models(SpreadInputs self)
#     cdef FuelModel get_fm_struct(SpreadInputs self, pyidx fm_number) noexcept
# cdef SpreadInputs make_SpreadInputs((float, float, float) cube_resolution, dict space_time_cubes)
# cdef CellInputs lookup_cell_inputs(SpreadInputs space_time_cubes, coord_tyx space_time_coordinate) noexcept
# cdef SpreadBehavior unburned_SpreadBehavior(vec_xy elevation_gradient, vec_xyz phi_gradient_xyz) noexcept
# cdef object encode_cell_index(pyidx y, pyidx x)
# cdef coord_yx decode_cell_index(object encoded_cell_index) noexcept
# cdef bint opposite_phi_signs(float[:,:] phi_matrix, pyidx y1, pyidx x1, pyidx y2, pyidx x2) noexcept
# cdef bint is_frontier_cell(
#     float[:,:] phi_matrix,
#     ISpaceTimeCube fuel_model_cube,
#     pyidx t,
#     pyidx y,
#     pyidx x,
#     ) noexcept
# cdef set identify_all_frontier_cells(
#     float[:,:] phi_matrix,
#     ISpaceTimeCube fuel_model_cube,
#     pyidx t,
#     pyidx rows,
#     pyidx cols,
#     )
# cdef NarrowBandTracker identify_tracked_cells(set frontier_cells, pyidx buffer_width, pyidx rows, pyidx cols)
# cdef void spot_from_burned_cell(
#     SpreadInputs space_time_cubes,
#     unsigned char[:,:] fire_type_matrix,
#     pyidx y,
#     pyidx x,
#     SpreadBehavior fire_behavior,
#     float time_of_arrival,
#     BufferedRandGen random_generator,
#     SpotConfig spot_config,
#     object spot_ignitions,
#     ) noexcept
# cdef float calc_phi_magnitude_xyz_2(vec_xy phi_gradient_xy, vec_xy elevation_gradient) noexcept
# cdef PartialedEllWavelet zero_partialed_wavelet() noexcept
# cdef PartialedEllWavelet prepare_partialed_wavelet(
#     vec_xyz heading_spread_vector,
#     float flanking_spread_rate,
#     float backing_spread_rate,
#     ) noexcept
# cdef PartialedEllWavelet wavelet_from_FireBehaviorMax(FireBehaviorMax fire_behavior_max) noexcept
# cdef float dphi_dt_from_partialed_wavelet(
#     PartialedEllWavelet wavelet,
#     vec_xy phi_gradient_xy,
#     float phi_magnitude_xyz_2,
#     ) noexcept
# cdef bint phi_aware_crowning_check(
#     float phi_magnitude_xyz_2,
#     float surface_dphi_dt,
#     float crowning_spread_rate,
#     ) noexcept
# cdef float dphi_dt_from_ellipses(EllipticalInfo ellipses, vec_xy phi_gradient_xy) noexcept
# cdef pyidx p_CellInputs = 17
# cdef class TrackedCellsArrays:
#     cdef pyidx _array_length
#     cdef pyidx n_tracked_cells
#     cdef float[:,:] float_inputs
#     cdef float[:,:] phi_values
#     cdef FireBehaviorMin* sfmin_arr
#     cdef EllipticalInfo* ell_info
#     cdef Pass1CellOutput* pass1outputs
#     cdef float[17] time_refreshed
#     cdef pyidx[17] t_refreshed
#     cdef void reset_size(TrackedCellsArrays self, pyidx n_tracked_cells)
# cdef void collect_phi_values(float[:,:] phi_matrix, TrackedCellsArrays tca) noexcept
# cdef int compare_cell_indexes(coord_yx c0, coord_yx c1) noexcept
# cdef void copy_tracked_cell_data(
#     pyidx i_old,
#     TrackedCellsArrays tca_old,
#     pyidx i_new,
#     TrackedCellsArrays tca_new,
#     ) noexcept
# cdef list[str] inputs_name_list()
# cdef class FireBehaviorSettings:
#     cdef float max_cells_per_timestep
#     cdef pyidx buffer_width
#     cdef bint use_wind_limit
#     cdef str surface_lw_ratio_model
#     cdef float crown_max_lw_ratio
#     cdef dict spot_config
#     cdef float[17] cube_refresh_rates
# cdef void load_float_inputs_for_cell(
#     SpreadInputs space_time_cubes,
#     coord_yx cell_index,
#     TrackedCellsArrays tca,
#     pyidx i,
#     ) noexcept
# cdef list[ISpaceTimeCube] list_float_input_cubes(SpreadInputs space_time_cubes)
# cdef dict default_cube_refresh_rates(float band_duration)
# cdef unsigned int recompute_level_for_input(pyidx input_k) noexcept
# cdef unsigned int refresh_inputs_if_needed(
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     TrackedCellsArrays tca,
#     float present_time,
#     ) noexcept
# cdef CellInputs load_saved_CellInputs(float[:,:] float_inputs, pyidx i) noexcept
# cdef FireBehaviorMin resolve_surface_no_wind_no_slope_behavior(CellInputs cell_inputs, FuelModel fuel_model) noexcept
# cdef FireBehaviorMax resolve_surface_max_behavior(
#     FireBehaviorSettings fb_opts,
#     CellInputs cell_inputs,
#     FuelModel fuel_model,
#     FireBehaviorMin surface_fire_min,
#     ) noexcept
# cdef FireBehaviorMax resolve_crown_max_behavior(
#     FireBehaviorSettings fb_opts,
#     CellInputs cell_inputs,
#     FuelModel fuel_model,
#     ) noexcept
# cdef float resolve_crowning_spread_rate(CellInputs cell_inputs, FireBehaviorMax surface_fire_max) noexcept
# cdef EllipticalInfo resolve_cell_elliptical_info(
#     FireBehaviorSettings fb_opts,
#     coord_yx cell_index,
#     CellInputs cell_inputs,
#     FuelModel fuel_model,
#     FireBehaviorMin surface_fire_min,
#     ) noexcept
# cdef void refresh_caches_from_inputs_if_needed(
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     TrackedCellsArrays tca,
#     float present_time,
#     ) noexcept
# cdef SpreadBehavior resolve_combined_spread_behavior(
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     coord_tyx space_time_coordinate,
#     vec_xy phi_gradient_xy,
#     ) noexcept
# cdef void load_tracked_cell_data(
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     coord_yx cell_index,
#     TrackedCellsArrays tca,
#     pyidx i,
#     ) noexcept
# cdef void sync_tracked_cells_arrays(
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     NarrowBandTracker tracked_cells,
#     TrackedCellsArrays tca_old,
#     TrackedCellsArrays tca_new,
#     ) noexcept
# cdef float runge_kutta_pass1(
#     float max_cells_per_timestep,
#     vec_xy spatial_resolution,
#     float max_timestep,
#     TrackedCellsArrays tca,
#     ) noexcept
# cdef void update_phi_star(TrackedCellsArrays tca, float dt, float[:,:] phi_star_matrix) noexcept
# cdef class BurnedCellInfo:
#     cdef coord_yx cell_index
#     cdef float time_of_arrival
#     cdef vec_xy phi_gradient_xy
#     cdef bint from_spotting
# cdef BurnedCellInfo new_BurnedCellInfo(
#     coord_yx cell_index,
#     float time_of_arrival,
#     vec_xy phi_gradient_xy,
#     bint from_spotting,
#     )
# cdef list[BurnedCellInfo] runge_kutta_pass2(
#     vec_xy spatial_resolution,
#     float start_time,
#     float dt,
#     TrackedCellsArrays tca,
#     float[:,:] phi_matrix,
#     )
# cdef void process_burned_cells(
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     dict output_matrices,
#     object spot_ignitions,
#     BufferedRandGen random_generator,
#     list[BurnedCellInfo] burned_cells,
#     ) noexcept
# cdef void reset_phi_star(
#     TrackedCellsArrays tca,
#     list[BurnedCellInfo] spot_ignited_cells,
#     float[:,:] phi_star_matrix,
#     float[:,:] phi_matrix,
#     ) noexcept
# cdef list[BurnedCellInfo] ignite_from_spotting(
#     object spot_ignitions,
#     dict output_matrices,
#     float stop_time,
#     )
# cdef void route_cell_to_diff(
#     set frontier_cells_old,
#     set frontier_additions,
#     set frontier_removals,
#     float[:,:] phi_matrix,
#     ISpaceTimeCube fuel_model_cube,
#     pyidx t,
#     pyidx y,
#     pyidx x,
#     ) noexcept
# cdef tuple[set, set] diff_frontier_cells(
#     set frontier_cells_old,
#     list[BurnedCellInfo] spread_ignited_cells,
#     list[BurnedCellInfo] spot_ignited_cells,
#     float[:,:] phi_matrix,
#     ISpaceTimeCube fuel_model_cube,
#     pyidx t,
#     )
# cdef set apply_frontier_diff(set frontier_cells_old, set frontier_additions, set frontier_removals)
# cdef NarrowBandTracker update_tracked_cells_with_frontier_diff(
#     NarrowBandTracker tracked_cells,
#     set frontier_cells_added,
#     set frontier_cells_dropped,
#     pyidx buffer_width,
#     )
# cdef dict spread_one_timestep(
#     dict sim_state,
#     SpreadInputs space_time_cubes,
#     FireBehaviorSettings fb_opts,
#     float max_timestep,
#     )
# cdef void check_space_time_cubes(dict space_time_cubes, object spot_config=?)
# cdef void check_output_matrices(dict output_matrices)
# cdef void check_dimensions_and_resolutions(
#     dict space_time_cubes,
#     dict output_matrices,
#     pyidx bands,
#     pyidx rows,
#     pyidx cols,
#     float band_duration,
#     float cell_height,
#     float cell_width,
#     )
# cdef void check_start_and_stop_times(float start_time, float max_stop_time, float cube_duration, object max_duration=?)
# cpdef dict[str, object] spread_fire_with_phi_field(
#     dict[str, ISpaceTimeCube] space_time_cubes,
#     dict[str, np.ndarray] output_matrices,
#     (float, float, float) cube_resolution,
#     float start_time,
#     object max_duration=?,
#     float max_cells_per_timestep=?,
#     pyidx buffer_width=?,
#     bint use_wind_limit=?,
#     str surface_lw_ratio_model=?,
#     float crown_max_lw_ratio=?,
#     dict[float, set] spot_ignitions=?,
#     object spot_config=?,
#     dict[str, float] cube_refresh_rates=?,
#     )
