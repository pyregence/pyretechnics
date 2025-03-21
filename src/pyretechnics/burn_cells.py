# [[file:../../org/pyretechnics.org::burn-cells-imports][burn-cells-imports]]
import cython
import cython as cy
import numpy as np
if cython.compiled:
    from cython.cimports.numpy import ndarray
    from cython.cimports.pyretechnics.cy_types import \
        pyidx, vec_xy, vec_xyz, coord_tyx, fclaarr, FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior
    from cython.cimports.pyretechnics.space_time_cube import ISpaceTimeCube, to_positive_index_range
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.vector_utils as vu
    import cython.cimports.pyretechnics.fuel_models as fm
    import cython.cimports.pyretechnics.surface_fire as sf
    import cython.cimports.pyretechnics.crown_fire as cf
else:
    from numpy import ndarray
    from pyretechnics.py_types import \
        pyidx, vec_xy, vec_xyz, coord_tyx, fclaarr, FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior
    from pyretechnics.space_time_cube import ISpaceTimeCube, to_positive_index_range
    import pyretechnics.conversion as conv
    import pyretechnics.vector_utils as vu
    import pyretechnics.fuel_models as fm
    import pyretechnics.surface_fire as sf
    import pyretechnics.crown_fire as cf
# burn-cells-imports ends here
# [[file:../../org/pyretechnics.org::burn-cells-as-head-fire][burn-cells-as-head-fire]]
@cy.cfunc
@cy.inline
def SpreadBehavior_to_dict(sb: SpreadBehavior) -> dict:
    return {
        "fire_type"         : sb.fire_type,
        "spread_rate"       : sb.spread_rate,
        "spread_direction"  : sb.spread_direction,
        "fireline_intensity": sb.fireline_intensity,
        "flame_length"      : sb.flame_length,
    }


@cy.ccall
def burn_cell_as_head_fire(space_time_cubes      : dict[str, ISpaceTimeCube],
                           space_time_coordinate : coord_tyx,
                           use_wind_limit        : cy.bint = True,
                           surface_lw_ratio_model: str = "behave",
                           crown_max_lw_ratio    : cy.float = 1e10) -> dict:
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_10hr       :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_100hr      :: kg moisture/kg ovendry weight
      - fuel_moisture_live_herbaceous :: kg moisture/kg ovendry weight
      - fuel_moisture_live_woody      :: kg moisture/kg ovendry weight
      - foliar_moisture               :: kg moisture/kg ovendry weight
      - fuel_spread_adjustment        :: float >= 0.0 (Optional: defaults to 1.0)
      - weather_spread_adjustment     :: float >= 0.0 (Optional: defaults to 1.0)
    - space_time_coordinate        :: (t,y,x)
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)

    return a dictionary with these fire behavior values for the space-time coordinate (t,y,x):
    - fire_type          :: 0 (unburned), 1 (surface), 2 (passive_crown), or 3 (active_crown)
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector on the slope-tangential plane
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    #================================================================================================
    # Destructure the space_time_coordinate
    #================================================================================================

    (t, y, x) = space_time_coordinate

    #================================================================================================
    # Unpack the space_time_cubes dictionary
    #================================================================================================

    # Topography, Fuel Model, and Vegetation
    slope              : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["slope"]).get(t, y, x)
    aspect             : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["aspect"]).get(t, y, x)
    fuel_model_number  : cy.int   = cy.cast(cy.int, cy.cast(ISpaceTimeCube, space_time_cubes["fuel_model"]).get(t, y, x))
    canopy_cover       : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_cover"]).get(t, y, x)
    canopy_height      : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_height"]).get(t, y, x)
    canopy_base_height : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_base_height"]).get(t, y, x)
    canopy_bulk_density: cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_bulk_density"]).get(t, y, x)

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m               : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["wind_speed_10m"]).get(t, y, x)
    upwind_direction             : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["upwind_direction"]).get(t, y, x)
    fuel_moisture_dead_1hr       : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_dead_1hr"]).get(t, y, x)
    fuel_moisture_dead_10hr      : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_dead_10hr"]).get(t, y, x)
    fuel_moisture_dead_100hr     : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_dead_100hr"]).get(t, y, x)
    fuel_moisture_live_herbaceous: cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_live_herbaceous"]).get(t, y, x)
    fuel_moisture_live_woody     : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_live_woody"]).get(t, y, x)
    foliar_moisture              : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["foliar_moisture"]).get(t, y, x)

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment   : cy.float = (cy.cast(ISpaceTimeCube, space_time_cubes["fuel_spread_adjustment"]).get(t, y, x)
                                           if "fuel_spread_adjustment" in space_time_cubes
                                           else 1.0)                                         # float >= 0.0
    weather_spread_adjustment: cy.float = (cy.cast(ISpaceTimeCube, space_time_cubes["weather_spread_adjustment"]).get(t, y, x)
                                           if "weather_spread_adjustment" in space_time_cubes
                                           else 1.0)                                         # float >= 0.0
    spread_rate_adjustment   : cy.float = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

    #================================================================================================
    # Check whether cell is burnable
    #================================================================================================

    fuel_model      : FuelModel
    maybe_fuel_model: FuelModel|None = fm.fuel_model_table.get(fuel_model_number)

    if maybe_fuel_model:
        fuel_model = maybe_fuel_model

    if (maybe_fuel_model is None or not(fuel_model.burnable)):
        # Encountered unknown or non-burnable fuel model

        #================================================================================================
        # Create a 3D unit vector pointing upslope on the slope-tangential plane
        #================================================================================================

        upslope_direction: cy.float = conv.opposite_direction(aspect)
        slope_vector_2d  : vec_xy   = conv.azimuthal_to_cartesian(slope, upslope_direction)
        slope_vector_3d  : vec_xyz  = vu.to_slope_plane(slope_vector_2d, slope_vector_2d)
        default_direction: vec_xyz  = (0.0, 1.0, 0.0) # default: North
        spread_direction : vec_xyz  = vu.as_unit_vector_3d(slope_vector_3d) if slope > 0.0 else default_direction

        #============================================================================================
        # Return zero surface fire behavior
        #============================================================================================

        return {
            "fire_type"         : 0, # unburned
            "spread_rate"       : 0.0,
            "spread_direction"  : spread_direction,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }

    else:
        # Encountered burnable fuel model

        #============================================================================================
        # Compute derived parameters
        #============================================================================================

        fuel_moisture               : fclaarr  = (fuel_moisture_dead_1hr,
                                                  fuel_moisture_dead_10hr,
                                                  fuel_moisture_dead_100hr,
                                                  0.0, # fuel_moisture_dead_herbaceous
                                                  fuel_moisture_live_herbaceous,
                                                  fuel_moisture_live_woody) # kg moisture/kg ovendry weight
        fuel_bed_depth              : cy.float = fuel_model.delta                      # ft
        heat_of_combustion          : cy.float = conv.Btu_lb_to_kJ_kg(fuel_model.h[0]) # kJ/kg
        estimated_fine_fuel_moisture: cy.float = fuel_moisture_dead_1hr                # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft: cy.float = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min: cy.float = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        midflame_wind_speed: cy.float = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,       # m/min
                                                                    fuel_bed_depth,              # ft
                                                                    conv.m_to_ft(canopy_height), # ft
                                                                    canopy_cover)                # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        moisturized_fuel_model: FuelModel = fm.moisturize(fuel_model, fuel_moisture)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min: FireBehaviorMin = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                                           spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max: FireBehaviorMax = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                                              midflame_wind_speed,
                                                                              upwind_direction,
                                                                              slope,
                                                                              aspect,
                                                                              use_wind_limit,
                                                                              surface_lw_ratio_model)

        # Simplify the surface fire behavior fields for future combination with the crown fire behavior values
        spread_direction       : vec_xyz        = surface_fire_max.max_spread_direction
        surface_fire_max_simple: SpreadBehavior = sf.calc_surface_fire_behavior_in_direction(surface_fire_max,
                                                                                             spread_direction)

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(surface_fire_max_simple.fireline_intensity,
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            # Calculate crown fire behavior in the direction of maximum spread
            crown_fire_max: FireBehaviorMax = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                                              canopy_bulk_density, heat_of_combustion,
                                                                              estimated_fine_fuel_moisture,
                                                                              wind_speed_10m, upwind_direction,
                                                                              slope, aspect, crown_max_lw_ratio)

            # Simplify the crown fire behavior fields for future combination with the surface fire behavior values
            spread_direction     : vec_xyz        = crown_fire_max.max_spread_direction
            crown_fire_max_simple: SpreadBehavior = cf.calc_crown_fire_behavior_in_direction(crown_fire_max,
                                                                                             spread_direction)

            #========================================================================================
            # Calculate combined fire behavior in the direction of maximum spread
            #========================================================================================

            combined_fire_max: SpreadBehavior = cf.calc_combined_fire_behavior(surface_fire_max_simple,
                                                                               crown_fire_max_simple)

            #========================================================================================
            # Return the combined fire behavior in the direction of maximum spread
            #========================================================================================

            return SpreadBehavior_to_dict(combined_fire_max)

        else:

            #========================================================================================
            # Return the surface fire behavior in the direction of maximum spread
            #========================================================================================

            return SpreadBehavior_to_dict(surface_fire_max_simple)


# TODO: Make a more efficient version that avoids space_time_cubes dictionary lookups for each cell
@cy.ccall
def burn_all_cells_as_head_fire(space_time_cubes      : dict[str, ISpaceTimeCube],
                                t                     : pyidx,
                                y_range               : tuple[pyidx, pyidx]|None = None,
                                x_range               : tuple[pyidx, pyidx]|None = None,
                                use_wind_limit        : cy.bint = True,
                                surface_lw_ratio_model: str = "behave",
                                crown_max_lw_ratio    : cy.float = 1e10) -> dict:
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_10hr       :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_100hr      :: kg moisture/kg ovendry weight
      - fuel_moisture_live_herbaceous :: kg moisture/kg ovendry weight
      - fuel_moisture_live_woody      :: kg moisture/kg ovendry weight
      - foliar_moisture               :: kg moisture/kg ovendry weight
      - fuel_spread_adjustment        :: float >= 0.0 (Optional: defaults to 1.0)
      - weather_spread_adjustment     :: float >= 0.0 (Optional: defaults to 1.0)
    - t                            :: temporal integer index into the ISpaceTimeCube objects
    - y_range                      :: (min_y, max_y) spatial integer index into the ISpaceTimeCube objects (Optional)
    - x_range                      :: (min_x, max_x) spatial integer index into the ISpaceTimeCube objects (Optional)
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)

    return a dictionary with these keys:
    - fire_type          :: 2D byte array (0=unburned, 1=surface, 2=passive_crown, 3=active_crown)
    - spread_rate        :: 2D float array (m/min)
    - spread_direction   :: 2D float array (degrees clockwise from North)
    - fireline_intensity :: 2D float array (kW/m)
    - flame_length       :: 2D float array (m)
    """
    slope_cube  : ISpaceTimeCube      = cy.cast(ISpaceTimeCube, space_time_cubes["slope"])
    bands       : pyidx               = slope_cube.shape[0]
    rows        : pyidx               = slope_cube.shape[1]
    cols        : pyidx               = slope_cube.shape[2]
    grid_shape  : tuple[int, int]     = (rows, cols)
    y_range_real: tuple[pyidx, pyidx] = to_positive_index_range(y_range, rows)
    x_range_real: tuple[pyidx, pyidx] = to_positive_index_range(x_range, cols)
    min_y       : pyidx               = y_range_real[0]
    max_y       : pyidx               = y_range_real[1]
    min_x       : pyidx               = x_range_real[0]
    max_x       : pyidx               = x_range_real[1]

    if not(0 <= t < bands):
        raise ValueError("The t value is out of range of the space_time_cubes.")

    if not(0 <= min_y < max_y <= rows):
        raise ValueError("The y_range values are out of range of the space_time_cubes.")

    if not(0 <= min_x < max_x <= rows):
        raise ValueError("The x_range values are out of range of the space_time_cubes.")

    fire_type_matrix          : ndarray         = np.zeros(grid_shape, dtype="uint8")
    spread_rate_matrix        : ndarray         = np.zeros(grid_shape, dtype="float32")
    spread_direction_matrix   : ndarray         = np.zeros(grid_shape, dtype="float32")
    fireline_intensity_matrix : ndarray         = np.zeros(grid_shape, dtype="float32")
    flame_length_matrix       : ndarray         = np.zeros(grid_shape, dtype="float32")
    fire_type_memview         : cy.uchar[:,::1] = fire_type_matrix
    spread_rate_memview       : cy.float[:,::1] = spread_rate_matrix
    spread_direction_memview  : cy.float[:,::1] = spread_direction_matrix
    fireline_intensity_memview: cy.float[:,::1] = fireline_intensity_matrix
    flame_length_memview      : cy.float[:,::1] = flame_length_matrix

    y                    : pyidx
    x                    : pyidx
    space_time_coordinate: coord_tyx
    for y in range(min_y, max_y):
        for x in range(min_x, max_x):
            space_time_coordinate          = (t, y, x)
            spread_behavior                = burn_cell_as_head_fire(space_time_cubes,
                                                                    space_time_coordinate,
                                                                    use_wind_limit,
                                                                    surface_lw_ratio_model,
                                                                    crown_max_lw_ratio)
            fire_type_memview[y,x]          = spread_behavior["fire_type"]
            spread_rate_memview[y,x]        = spread_behavior["spread_rate"]
            spread_direction_memview[y,x]   = vu.spread_direction_vector_to_angle(spread_behavior["spread_direction"])
            fireline_intensity_memview[y,x] = spread_behavior["fireline_intensity"]
            flame_length_memview[y,x]       = spread_behavior["flame_length"]

    return {
        "fire_type"         : fire_type_matrix,
        "spread_rate"       : spread_rate_matrix,
        "spread_direction"  : spread_direction_matrix,
        "fireline_intensity": fireline_intensity_matrix,
        "flame_length"      : flame_length_matrix,
    }
# burn-cells-as-head-fire ends here
# [[file:../../org/pyretechnics.org::burn-cells-toward-azimuth][burn-cells-toward-azimuth]]
@cy.ccall
def burn_cell_toward_azimuth(space_time_cubes      : dict[str, ISpaceTimeCube],
                             space_time_coordinate : coord_tyx,
                             azimuth               : cy.float,
                             use_wind_limit        : cy.bint = True,
                             surface_lw_ratio_model: str = "behave",
                             crown_max_lw_ratio    : cy.float = 1e10) -> dict:
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_10hr       :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_100hr      :: kg moisture/kg ovendry weight
      - fuel_moisture_live_herbaceous :: kg moisture/kg ovendry weight
      - fuel_moisture_live_woody      :: kg moisture/kg ovendry weight
      - foliar_moisture               :: kg moisture/kg ovendry weight
      - fuel_spread_adjustment        :: float >= 0.0 (Optional: defaults to 1.0)
      - weather_spread_adjustment     :: float >= 0.0 (Optional: defaults to 1.0)
    - space_time_coordinate        :: (t,y,x)
    - azimuth                      :: degrees clockwise from North on the horizontal plane
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)

    return a dictionary with these fire behavior values for the space-time coordinate (t,y,x):
    - fire_type          :: 0 (unburned), 1 (surface), 2 (passive_crown), or 3 (active_crown)
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector on the slope-tangential plane
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    #================================================================================================
    # Destructure the space_time_coordinate
    #================================================================================================

    (t, y, x) = space_time_coordinate

    #================================================================================================
    # Unpack the space_time_cubes dictionary
    #================================================================================================

    # Topography, Fuel Model, and Vegetation
    slope              : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["slope"]).get(t, y, x)
    aspect             : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["aspect"]).get(t, y, x)
    fuel_model_number  : cy.int   = cy.cast(cy.int, cy.cast(ISpaceTimeCube, space_time_cubes["fuel_model"]).get(t, y, x))
    canopy_cover       : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_cover"]).get(t, y, x)
    canopy_height      : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_height"]).get(t, y, x)
    canopy_base_height : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_base_height"]).get(t, y, x)
    canopy_bulk_density: cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["canopy_bulk_density"]).get(t, y, x)

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m               : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["wind_speed_10m"]).get(t, y, x)
    upwind_direction             : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["upwind_direction"]).get(t, y, x)
    fuel_moisture_dead_1hr       : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_dead_1hr"]).get(t, y, x)
    fuel_moisture_dead_10hr      : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_dead_10hr"]).get(t, y, x)
    fuel_moisture_dead_100hr     : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_dead_100hr"]).get(t, y, x)
    fuel_moisture_live_herbaceous: cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_live_herbaceous"]).get(t, y, x)
    fuel_moisture_live_woody     : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["fuel_moisture_live_woody"]).get(t, y, x)
    foliar_moisture              : cy.float = cy.cast(ISpaceTimeCube, space_time_cubes["foliar_moisture"]).get(t, y, x)

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment   : cy.float = (cy.cast(ISpaceTimeCube, space_time_cubes["fuel_spread_adjustment"]).get(t, y, x)
                                           if "fuel_spread_adjustment" in space_time_cubes
                                           else 1.0)                                         # float >= 0.0
    weather_spread_adjustment: cy.float = (cy.cast(ISpaceTimeCube, space_time_cubes["weather_spread_adjustment"]).get(t, y, x)
                                           if "weather_spread_adjustment" in space_time_cubes
                                           else 1.0)                                         # float >= 0.0
    spread_rate_adjustment   : cy.float = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

    #================================================================================================
    # Project a 2D unit vector pointing toward the azimuth onto the slope-tangential plane
    #================================================================================================

    upslope_direction: cy.float = conv.opposite_direction(aspect)
    slope_vector_2d  : vec_xy   = conv.azimuthal_to_cartesian(slope, upslope_direction)
    azimuth_vector_2d: vec_xy   = conv.azimuthal_to_cartesian(1.0, azimuth)
    spread_direction : vec_xyz  = vu.as_unit_vector_3d(vu.to_slope_plane(azimuth_vector_2d, slope_vector_2d))

    #================================================================================================
    # Check whether cell is burnable
    #================================================================================================

    fuel_model      : FuelModel
    maybe_fuel_model: FuelModel|None = fm.fuel_model_table.get(fuel_model_number)

    if maybe_fuel_model:
        fuel_model = maybe_fuel_model

    if (maybe_fuel_model is None or not(fuel_model.burnable)):
        # Encountered unknown or non-burnable fuel model

        #============================================================================================
        # Return zero surface fire behavior in the direction of the azimuth vector
        #============================================================================================

        return {
            "fire_type"         : 0, # unburned
            "spread_rate"       : 0.0,
            "spread_direction"  : spread_direction,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }

    else:
        # Encountered burnable fuel model

        #============================================================================================
        # Compute derived parameters
        #============================================================================================

        fuel_moisture               : fclaarr  = (fuel_moisture_dead_1hr,
                                                  fuel_moisture_dead_10hr,
                                                  fuel_moisture_dead_100hr,
                                                  0.0, # fuel_moisture_dead_herbaceous
                                                  fuel_moisture_live_herbaceous,
                                                  fuel_moisture_live_woody) # kg moisture/kg ovendry weight
        fuel_bed_depth              : cy.float = fuel_model.delta                      # ft
        heat_of_combustion          : cy.float = conv.Btu_lb_to_kJ_kg(fuel_model.h[0]) # kJ/kg
        estimated_fine_fuel_moisture: cy.float = fuel_moisture_dead_1hr                # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft: cy.float = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min: cy.float = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        midflame_wind_speed: cy.float = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,       # m/min
                                                                    fuel_bed_depth,              # ft
                                                                    conv.m_to_ft(canopy_height), # ft
                                                                    canopy_cover)                # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        moisturized_fuel_model: FuelModel = fm.moisturize(fuel_model, fuel_moisture)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min: FireBehaviorMin = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                                           spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max: FireBehaviorMax = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                                              midflame_wind_speed,
                                                                              upwind_direction,
                                                                              slope,
                                                                              aspect,
                                                                              use_wind_limit,
                                                                              surface_lw_ratio_model)

        #============================================================================================
        # Calculate surface fire behavior in the direction of the azimuth vector
        #============================================================================================

        surface_fire_azimuth: SpreadBehavior = sf.calc_surface_fire_behavior_in_direction(surface_fire_max,
                                                                                          spread_direction)

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(surface_fire_azimuth.fireline_intensity,
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            crown_fire_max: FireBehaviorMax = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                                              canopy_bulk_density, heat_of_combustion,
                                                                              estimated_fine_fuel_moisture,
                                                                              wind_speed_10m, upwind_direction,
                                                                              slope, aspect, crown_max_lw_ratio)

            #========================================================================================
            # Calculate crown fire behavior in the direction of the azimuth vector
            #========================================================================================

            crown_fire_azimuth: SpreadBehavior = cf.calc_crown_fire_behavior_in_direction(crown_fire_max,
                                                                                          spread_direction)

            #========================================================================================
            # Calculate combined fire behavior in the direction of the azimuth vector
            #========================================================================================

            combined_fire_azimuth: SpreadBehavior = cf.calc_combined_fire_behavior(surface_fire_azimuth,
                                                                                   crown_fire_azimuth)

            #========================================================================================
            # Return the combined fire behavior in the direction of the azimuth vector
            #========================================================================================

            return SpreadBehavior_to_dict(combined_fire_azimuth)

        else:

            #========================================================================================
            # Return the surface fire behavior in the direction of the azimuth vector
            #========================================================================================

            return SpreadBehavior_to_dict(surface_fire_azimuth)


# TODO: Make a more efficient version that avoids space_time_cubes dictionary lookups for each cell
@cy.ccall
def burn_all_cells_toward_azimuth(space_time_cubes      : dict[str, ISpaceTimeCube],
                                  azimuth               : cy.float,
                                  t                     : pyidx,
                                  y_range               : tuple[pyidx, pyidx]|None = None,
                                  x_range               : tuple[pyidx, pyidx]|None = None,
                                  use_wind_limit        : cy.bint = True,
                                  surface_lw_ratio_model: str = "behave",
                                  crown_max_lw_ratio    : cy.float = 1e10) -> dict:
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_10hr       :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_100hr      :: kg moisture/kg ovendry weight
      - fuel_moisture_live_herbaceous :: kg moisture/kg ovendry weight
      - fuel_moisture_live_woody      :: kg moisture/kg ovendry weight
      - foliar_moisture               :: kg moisture/kg ovendry weight
      - fuel_spread_adjustment        :: float >= 0.0 (Optional: defaults to 1.0)
      - weather_spread_adjustment     :: float >= 0.0 (Optional: defaults to 1.0)
    - azimuth                      :: degrees clockwise from North on the horizontal plane
    - t                            :: temporal integer index into the ISpaceTimeCube objects
    - y_range                      :: (min_y, max_y) spatial integer index into the ISpaceTimeCube objects (Optional)
    - x_range                      :: (min_x, max_x) spatial integer index into the ISpaceTimeCube objects (Optional)
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)

    return a dictionary with these keys:
    - fire_type          :: 2D byte array (0=unburned, 1=surface, 2=passive_crown, 3=active_crown)
    - spread_rate        :: 2D float array (m/min)
    - spread_direction   :: 2D float array (degrees clockwise from North)
    - fireline_intensity :: 2D float array (kW/m)
    - flame_length       :: 2D float array (m)
    """
    slope_cube  : ISpaceTimeCube      = cy.cast(ISpaceTimeCube, space_time_cubes["slope"])
    bands       : pyidx               = slope_cube.shape[0]
    rows        : pyidx               = slope_cube.shape[1]
    cols        : pyidx               = slope_cube.shape[2]
    grid_shape  : tuple[int, int]     = (rows, cols)
    y_range_real: tuple[pyidx, pyidx] = to_positive_index_range(y_range, rows)
    x_range_real: tuple[pyidx, pyidx] = to_positive_index_range(x_range, cols)
    min_y       : pyidx               = y_range_real[0]
    max_y       : pyidx               = y_range_real[1]
    min_x       : pyidx               = x_range_real[0]
    max_x       : pyidx               = x_range_real[1]

    if not(0 <= t < bands):
        raise ValueError("The t value is out of range of the space_time_cubes.")

    if not(0 <= min_y < max_y <= rows):
        raise ValueError("The y_range values are out of range of the space_time_cubes.")

    if not(0 <= min_x < max_x <= rows):
        raise ValueError("The x_range values are out of range of the space_time_cubes.")

    fire_type_matrix          : ndarray         = np.zeros(grid_shape, dtype="uint8")
    spread_rate_matrix        : ndarray         = np.zeros(grid_shape, dtype="float32")
    spread_direction_matrix   : ndarray         = np.zeros(grid_shape, dtype="float32")
    fireline_intensity_matrix : ndarray         = np.zeros(grid_shape, dtype="float32")
    flame_length_matrix       : ndarray         = np.zeros(grid_shape, dtype="float32")
    fire_type_memview         : cy.uchar[:,::1] = fire_type_matrix
    spread_rate_memview       : cy.float[:,::1] = spread_rate_matrix
    spread_direction_memview  : cy.float[:,::1] = spread_direction_matrix
    fireline_intensity_memview: cy.float[:,::1] = fireline_intensity_matrix
    flame_length_memview      : cy.float[:,::1] = flame_length_matrix

    y                    : pyidx
    x                    : pyidx
    space_time_coordinate: coord_tyx
    for y in range(min_y, max_y):
        for x in range(min_x, max_x):
            space_time_coordinate          = (t, y, x)
            spread_behavior                = burn_cell_toward_azimuth(space_time_cubes,
                                                                      space_time_coordinate,
                                                                      azimuth,
                                                                      use_wind_limit,
                                                                      surface_lw_ratio_model,
                                                                      crown_max_lw_ratio)
            fire_type_memview[y,x]          = spread_behavior["fire_type"]
            spread_rate_memview[y,x]        = spread_behavior["spread_rate"]
            spread_direction_memview[y,x]   = vu.spread_direction_vector_to_angle(spread_behavior["spread_direction"])
            fireline_intensity_memview[y,x] = spread_behavior["fireline_intensity"]
            flame_length_memview[y,x]       = spread_behavior["flame_length"]

    return {
        "fire_type"         : fire_type_matrix,
        "spread_rate"       : spread_rate_matrix,
        "spread_direction"  : spread_direction_matrix,
        "fireline_intensity": fireline_intensity_matrix,
        "flame_length"      : flame_length_matrix,
    }
# burn-cells-toward-azimuth ends here
