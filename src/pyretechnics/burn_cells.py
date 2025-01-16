# [[file:../../org/pyretechnics.org::burn-cell-as-head-fire][burn-cell-as-head-fire]]
# cython: profile=False
import cython
if cython.compiled:
    from cython.cimports.pyretechnics.cy_types import vec_xy, vec_xyz
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.vector_utils as vu
    import cython.cimports.pyretechnics.fuel_models as fm
    import cython.cimports.pyretechnics.surface_fire as sf
    import cython.cimports.pyretechnics.crown_fire as cf
else:
    from pyretechnics.py_types import vec_xy, vec_xyz
    import pyretechnics.conversion as conv
    import pyretechnics.vector_utils as vu
    import pyretechnics.fuel_models as fm
    import pyretechnics.surface_fire as sf
    import pyretechnics.crown_fire as cf


import cython as cy


# TODO: Create a version of this function that runs efficiently over a space_time_region
def burn_cell_as_head_fire(space_time_cubes, space_time_coordinate, use_wind_limit=True,
                           surface_lw_ratio_model="rothermel", crown_max_lw_ratio=None):
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
    - fire_type          :: "unburned", "surface", "passive_crown", or "active_crown"
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
    slope               = space_time_cubes["slope"].get(t,y,x)               # rise/run
    aspect              = space_time_cubes["aspect"].get(t,y,x)              # degrees clockwise from North
    fuel_model_number   = space_time_cubes["fuel_model"].get(t,y,x)          # integer index in fm.fuel_model_table
    canopy_cover        = space_time_cubes["canopy_cover"].get(t,y,x)        # 0-1
    canopy_height       = space_time_cubes["canopy_height"].get(t,y,x)       # m
    canopy_base_height  = space_time_cubes["canopy_base_height"].get(t,y,x)  # m
    canopy_bulk_density = space_time_cubes["canopy_bulk_density"].get(t,y,x) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                = space_time_cubes["wind_speed_10m"].get(t,y,x)                # km/hr
    upwind_direction              = space_time_cubes["upwind_direction"].get(t,y,x)              # degrees clockwise from North
    fuel_moisture_dead_1hr        = space_time_cubes["fuel_moisture_dead_1hr"].get(t,y,x)        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       = space_time_cubes["fuel_moisture_dead_10hr"].get(t,y,x)       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      = space_time_cubes["fuel_moisture_dead_100hr"].get(t,y,x)      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous = space_time_cubes["fuel_moisture_live_herbaceous"].get(t,y,x) # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      = space_time_cubes["fuel_moisture_live_woody"].get(t,y,x)      # kg moisture/kg ovendry weight
    foliar_moisture               = space_time_cubes["foliar_moisture"].get(t,y,x)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    = (space_time_cubes["fuel_spread_adjustment"].get(t,y,x)
                                 if "fuel_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    weather_spread_adjustment = (space_time_cubes["weather_spread_adjustment"].get(t,y,x)
                                 if "weather_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    spread_rate_adjustment    = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

    #================================================================================================
    # Check whether cell is burnable
    #================================================================================================

    fuel_model = fm.fuel_model_table.get(fuel_model_number)

    if not (fuel_model and fuel_model["burnable"]):
        # Encountered unknown or non-burnable fuel model

        #================================================================================================
        # Create a 3D unit vector pointing upslope on the slope-tangential plane
        #================================================================================================

        upslope_direction: cy.float = conv.opposite_direction(aspect)
        slope_vector_2d  : vec_xy   = conv.azimuthal_to_cartesian(slope, upslope_direction)
        slope_vector_3d  : vec_xyz  = vu.to_slope_plane(slope_vector_2d, slope_vector_2d)
        spread_direction : vec_xyz  = vu.as_unit_vector_3d(slope_vector_3d) if slope > 0.0 else (0.0,1.0,0.0) # default: North

        #============================================================================================
        # Return zero surface fire behavior
        #============================================================================================

        return {
            "fire_type"         : "unburned",
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

        fuel_moisture                = [fuel_moisture_dead_1hr,
                                        fuel_moisture_dead_10hr,
                                        fuel_moisture_dead_100hr,
                                        0.0, # fuel_moisture_dead_herbaceous
                                        fuel_moisture_live_herbaceous,
                                        fuel_moisture_live_woody]               # kg moisture/kg ovendry weight
        fuel_bed_depth               = fuel_model["delta"]                      # ft
        heat_of_combustion           = conv.Btu_lb_to_kJ_kg(fuel_model["h"][0]) # kJ/kg
        estimated_fine_fuel_moisture = fuel_moisture_dead_1hr                   # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        midflame_wind_speed = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,       # m/min
                                                          fuel_bed_depth,              # ft
                                                          conv.m_to_ft(canopy_height), # ft
                                                          canopy_cover)                # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        moisturized_fuel_model = fm.moisturize(fuel_model, fuel_moisture)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                          spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                             midflame_wind_speed,
                                                             upwind_direction,
                                                             slope,
                                                             aspect,
                                                             use_wind_limit,
                                                             surface_lw_ratio_model)

        # Simplify the surface fire behavior fields for future comparison/combination with the crown fire behavior values
        surface_fire_max_simple = sf.calc_surface_fire_behavior_in_direction(surface_fire_max,
                                                                             surface_fire_max["max_spread_direction"])

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(surface_fire_max_simple["fireline_intensity"],
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            # Calculate crown fire behavior in the direction of maximum spread
            crown_fire_max = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                             canopy_bulk_density, heat_of_combustion,
                                                             estimated_fine_fuel_moisture,
                                                             wind_speed_10m, upwind_direction,
                                                             slope, aspect, crown_max_lw_ratio)

            # Simplify the crown fire behavior fields for future comparison/combination with the surface fire behavior values
            crown_fire_max_simple = cf.calc_crown_fire_behavior_in_direction(crown_fire_max,
                                                                             crown_fire_max.max_spread_direction)

            #========================================================================================
            # Calculate combined fire behavior in the direction of maximum spread
            #========================================================================================

            combined_fire_max = cf.calc_combined_fire_behavior(surface_fire_max_simple, crown_fire_max_simple)

            #========================================================================================
            # Return the combined fire behavior in the direction of maximum spread
            #========================================================================================

            return combined_fire_max

        else:

            #========================================================================================
            # Return the surface fire behavior in the direction of maximum spread
            #========================================================================================

            return surface_fire_max_simple
# burn-cell-as-head-fire ends here
# [[file:../../org/pyretechnics.org::burn-cell-toward-azimuth][burn-cell-toward-azimuth]]
# TODO: Create a version of this function that runs efficiently over a space_time_region
def burn_cell_toward_azimuth(space_time_cubes, space_time_coordinate, azimuth, use_wind_limit=True,
                             surface_lw_ratio_model="rothermel", crown_max_lw_ratio=None):
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
    - fire_type          :: "unburned", "surface", "passive_crown", or "active_crown"
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
    slope               = space_time_cubes["slope"].get(t,y,x)               # rise/run
    aspect              = space_time_cubes["aspect"].get(t,y,x)              # degrees clockwise from North
    fuel_model_number   = space_time_cubes["fuel_model"].get(t,y,x)          # integer index in fm.fuel_model_table
    canopy_cover        = space_time_cubes["canopy_cover"].get(t,y,x)        # 0-1
    canopy_height       = space_time_cubes["canopy_height"].get(t,y,x)       # m
    canopy_base_height  = space_time_cubes["canopy_base_height"].get(t,y,x)  # m
    canopy_bulk_density = space_time_cubes["canopy_bulk_density"].get(t,y,x) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                = space_time_cubes["wind_speed_10m"].get(t,y,x)                # km/hr
    upwind_direction              = space_time_cubes["upwind_direction"].get(t,y,x)              # degrees clockwise from North
    fuel_moisture_dead_1hr        = space_time_cubes["fuel_moisture_dead_1hr"].get(t,y,x)        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       = space_time_cubes["fuel_moisture_dead_10hr"].get(t,y,x)       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      = space_time_cubes["fuel_moisture_dead_100hr"].get(t,y,x)      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous = space_time_cubes["fuel_moisture_live_herbaceous"].get(t,y,x) # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      = space_time_cubes["fuel_moisture_live_woody"].get(t,y,x)      # kg moisture/kg ovendry weight
    foliar_moisture               = space_time_cubes["foliar_moisture"].get(t,y,x)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    = (space_time_cubes["fuel_spread_adjustment"].get(t,y,x)
                                 if "fuel_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    weather_spread_adjustment = (space_time_cubes["weather_spread_adjustment"].get(t,y,x)
                                 if "weather_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    spread_rate_adjustment    = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

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

    fuel_model = fm.fuel_model_table.get(fuel_model_number)

    if not (fuel_model and fuel_model["burnable"]):
        # Encountered unknown or non-burnable fuel model

        #============================================================================================
        # Return zero surface fire behavior in the direction of the azimuth vector
        #============================================================================================

        return {
            "fire_type"         : "unburned",
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

        fuel_moisture                = [fuel_moisture_dead_1hr,
                                        fuel_moisture_dead_10hr,
                                        fuel_moisture_dead_100hr,
                                        0.0, # fuel_moisture_dead_herbaceous
                                        fuel_moisture_live_herbaceous,
                                        fuel_moisture_live_woody]               # kg moisture/kg ovendry weight
        fuel_bed_depth               = fuel_model["delta"]                      # ft
        heat_of_combustion           = conv.Btu_lb_to_kJ_kg(fuel_model["h"][0]) # kJ/kg
        estimated_fine_fuel_moisture = fuel_moisture_dead_1hr                   # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        midflame_wind_speed = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,       # m/min
                                                          fuel_bed_depth,              # ft
                                                          conv.m_to_ft(canopy_height), # ft
                                                          canopy_cover)                # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        moisturized_fuel_model = fm.moisturize(fuel_model, fuel_moisture)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                          spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                             midflame_wind_speed,
                                                             upwind_direction,
                                                             slope,
                                                             aspect,
                                                             use_wind_limit,
                                                             surface_lw_ratio_model)

        #============================================================================================
        # Calculate surface fire behavior in the direction of the azimuth vector
        #============================================================================================

        surface_fire_azimuth = sf.calc_surface_fire_behavior_in_direction(surface_fire_max, spread_direction)

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(surface_fire_azimuth["fireline_intensity"],
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            crown_fire_max = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                             canopy_bulk_density, heat_of_combustion,
                                                             estimated_fine_fuel_moisture,
                                                             wind_speed_10m, upwind_direction,
                                                             slope, aspect, crown_max_lw_ratio)

            #========================================================================================
            # Calculate crown fire behavior in the direction of the azimuth vector
            #========================================================================================

            crown_fire_azimuth = cf.calc_crown_fire_behavior_in_direction(crown_fire_max, spread_direction)

            #========================================================================================
            # Calculate combined fire behavior in the direction of the azimuth vector
            #========================================================================================

            combined_fire_azimuth = cf.calc_combined_fire_behavior(surface_fire_azimuth, crown_fire_azimuth)

            #========================================================================================
            # Return the combined fire behavior in the direction of the azimuth vector
            #========================================================================================

            return combined_fire_azimuth

        else:

            #========================================================================================
            # Return the surface fire behavior in the direction of the azimuth vector
            #========================================================================================

            return surface_fire_azimuth
# burn-cell-toward-azimuth ends here
