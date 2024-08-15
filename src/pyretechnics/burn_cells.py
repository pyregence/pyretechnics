# [[file:../../org/pyretechnics.org::burn-cells][burn-cells]]
import pyretechnics.conversion as conv
import pyretechnics.fuel_models as fm
import pyretechnics.surface_fire as sf
import pyretechnics.crown_fire as cf


def one_everywhere(t, y, x):
    """
    Return 1.0 for all t, y, x values.
    """
    return 1.0


# NOTE: Name change (no more !)
def compute_max_in_situ_values(inputs, t, y, x):
    """
    Returns the following fire behavior values for the space-time cell at location (t,y,x):
    - max_spread_rate       : m/min
    - max_spread_direction  : deg
    - max_fireline_intensity: kW/m
    - max_flame_length      : m
    - fire_type             : 0=unburned, 1=surface, 2=passive crown, 3=active crown
    - eccentricity          : unitless (0: circular fire, >0: elliptical fire)
    """
    # Topography, Fuel Model, and Vegetation
    slope                         = inputs["slope"](t,y,x)
    aspect                        = inputs["aspect"](t,y,x)
    fuel_model_number             = inputs["fuel_model"](t,y,x)
    canopy_cover                  = inputs["canopy_cover"](t,y,x)
    canopy_height                 = inputs["canopy_height"](t,y,x)
    canopy_base_height            = inputs["canopy_base_height"](t,y,x)
    canopy_bulk_density           = inputs["canopy_bulk_density"](t,y,x)
    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m_x              = inputs["wind_speed_10m_x"](t,y,x)
    wind_speed_10m_y              = inputs["wind_speed_10m_y"](t,y,x)
    fuel_moisture_dead_1hr        = inputs["fuel_moisture_dead_1hr"](t,y,x)
    fuel_moisture_dead_10hr       = inputs["fuel_moisture_dead_10hr"](t,y,x)
    fuel_moisture_dead_100hr      = inputs["fuel_moisture_dead_100hr"](t,y,x)
    fuel_moisture_live_herbaceous = inputs["fuel_moisture_live_herbaceous"](t,y,x)
    fuel_moisture_live_woody      = inputs["fuel_moisture_live_woody"](t,y,x)
    foliar_moisture               = inputs["foliar_moisture"](t,y,x)
    # Spread Rate Adjustments
    fuel_spread_adjustment        = inputs.get("fuel_spread_adjustment"   , one_everywhere)(t,y,x)
    weather_spread_adjustment     = inputs.get("weather_spread_adjustment", one_everywhere)(t,y,x)
    spread_rate_adjustment        = fuel_spread_adjustment * weather_spread_adjustment
    # Check Whether Cell is Burnable
    fuel_model = fm.fuel_model_table.get(fuel_model_number)
    if not (fuel_model and fuel_model["burnable"]):
        return {
            "max_spread_rate"       : 0.0, # m/min
            "max_spread_direction"  : 0.0, # deg
            "max_fireline_intensity": 0.0, # kW/m
            "max_flame_length"      : 0.0, # m
            "fire_type"             : 0,   # 0=unburned
            "eccentricity"          : 0.0, # unitless
        }
    else:
        # Moisturized Fuel Model
        fuel_moisture           = [fuel_moisture_dead_1hr,
                                   fuel_moisture_dead_10hr,
                                   fuel_moisture_dead_100hr,
                                   0.0, # fuel_moisture_dead_herbaceous
                                   fuel_moisture_live_herbaceous,
                                   fuel_moisture_live_woody]
        moisturized_fuel_model  = fm.moisturize(fuel_model, fuel_moisture)
        # Baseline Surface Fire Behavior
        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        surface_fire_min        = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                                 spread_rate_adjustment)
        # Midflame Wind Speed
        (wind_speed_10m,
         upwind_direction)      = conv.cartesian_to_azimuthal(wind_speed_10m_x, wind_speed_10m_y) # (m/min, deg)
        wind_speed_20ft         = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # m/min
        midflame_wind_speed     = sf.calc_midflame_wind_speed(conv.m_to_ft(wind_speed_20ft),
                                                              fuel_model["delta"],
                                                              conv.m_to_ft(canopy_height),
                                                              canopy_cover) # ft/min
        # Max Surface Fire Behavior
        surface_fire_max        = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                                    midflame_wind_speed,
                                                                    upwind_direction,
                                                                    slope,
                                                                    aspect)
        max_surface_spread_rate      = surface_fire_max["max_spread_rate"] # ft/min
        max_spread_direction         = surface_fire_max["max_spread_direction"] # deg
        surface_eccentricity         = surface_fire_max["eccentricity"] # unitless
        max_surface_intensity        = surface_fire_max["max_fireline_intensity"] # Btu/ft/s
        max_surface_intensity_metric = conv.Btu_ft_s_to_kW_m(max_surface_intensity) # kW/m
        #=======================================================================================
        # TODO: The calculations to determine the fireline_normal_spread_rate have been elided.
        #       Consider ending this function here and making another function to compute the
        #       surface/crown values based on the provided perimeter spread direction.
        #=======================================================================================
        # Check for Crown Fire Initiation
        if cf.van_wagner_crown_fire_initiation(max_surface_intensity_metric,
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):
            # Max Crown Spread Rate, Fire Type, and Crown Eccentricity
            crown_spread_info     = cf.cruz_crown_fire_spread_info(conv.m_min_to_km_hr(wind_speed_10m),
                                                                   canopy_bulk_density,
                                                                   fuel_moisture_dead_1hr)
            fire_type             = crown_spread_info["fire_type"] # "passive" or "active"
            max_crown_spread_rate = conv.m_to_ft(crown_spread_info["spread_rate"]) # ft/min
            crown_eccentricity    = (
                surface_eccentricity
                if max_surface_spread_rate > max_crown_spread_rate
                else cf.crown_fire_eccentricity(cf.crown_length_to_width_ratio(conv.m_min_to_mph(wind_speed_20ft)))
            ) # unitless
            #=======================================================================================
            # NOTE: The calculations to determine the fireline_normal_spread_rate have been elided.
            #=======================================================================================
            # Max Crown Intensity
            max_crown_intensity   = cf.calc_crown_fireline_intensity(conv.ft_to_m(max_crown_spread_rate),
                                                                     canopy_bulk_density,
                                                                     canopy_height,
                                                                     canopy_base_height,
                                                                     conv.Btu_lb_to_kJ_kg(fuel_model["h"][0])) # kW/m
            # Max Combined Spread Rate, Intensity, and Flame Length
            max_spread_rate       = conv.ft_to_m(max(max_surface_spread_rate, max_crown_spread_rate)) # m/min
            max_intensity         = max_surface_intensity_metric + max_crown_intensity # kW/m
            max_flame_length      = conv.ft_to_m(sf.calc_flame_length(conv.kW_m_to_Btu_ft_s(max_intensity))) # m
            # Return Fire Behavior Values
            return {
                "max_spread_rate"       : max_spread_rate,      # m/min
                "max_spread_direction"  : max_spread_direction, # deg
                "max_fireline_intensity": max_intensity,        # kW/m
                "max_flame_length"      : max_flame_length,     # m
                "fire_type"             : fire_type,            # "passive" or "active"
                "eccentricity"          : crown_eccentricity,   # unitless
            }
        else:
            max_surface_flame_length = conv.ft_to_m(sf.calc_flame_length(max_surface_intensity)) # m
            return {
                "max_spread_rate"       : conv.ft_to_m(max_surface_spread_rate), # m/min
                "max_spread_direction"  : max_spread_direction,                  # deg
                "max_fireline_intensity": max_surface_intensity_metric,          # kW/m
                "max_flame_length"      : max_surface_flame_length,              # m
                "fire_type"             : "surface",                             # "surface"
                "eccentricity"          : surface_eccentricity,                  # unitless
            }
# burn-cells ends here
