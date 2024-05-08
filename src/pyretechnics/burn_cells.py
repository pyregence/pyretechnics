# [[file:../../org/pyretechnics.org::burn-cells][burn-cells]]
from math import sqrt, atan2, degrees
from pyretechnics.conversion import wind_speed_10m_to_wind_speed_20ft, m_to_ft, Btu_ft_s_to_kW_m, kW_m_to_Btu_ft_s, m_min_to_km_hr, m_min_to_mph, ft_to_m, Btu_lb_to_kJ_kg
from pyretechnics.fuel_models import fuel_models_precomputed, moisturize, is_burnable_fuel_model_number
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
    - max_spread_rate        : m/min
    - max_spread_direction   : deg
    - max_fire_line_intensity: kW/m
    - max_flame_length       : m
    - fire_type              : 0=unburned, 1=surface, 2=passive crown, 3=active crown
    - eccentricity           : unitless (0: circular fire, >0: elliptical fire)
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
    # Check Whether Cell is Burnable
    if not is_burnable_fuel_model_number(fuel_model_number):
        return {
            "max_spread_rate"        : 0.0, # m/min
            "max_spread_direction"   : 0.0, # deg
            "max_fire_line_intensity": 0.0, # kW/m
            "max_flame_length"       : 0.0, # m
            "fire_type"              : 0,   # 0=unburned
            "eccentricity"           : 0.0, # unitless
        }
    else:
        # Moisturized Fuel Model
        fuel_model              = fuel_models_precomputed[fuel_model_number]
        fuel_moisture           = [fuel_moisture_dead_1hr,
                                   fuel_moisture_dead_10hr,
                                   fuel_moisture_dead_100hr,
                                   0.0, # fuel_moisture_dead_herbaceous
                                   fuel_moisture_live_herbaceous,
                                   fuel_moisture_live_woody]
        moisturized_fuel_model  = moisturize(fuel_model, fuel_moisture)
        # Baseline Surface Spread Rate, Residence Time, and Reaction Intensity
        # TODO: Memoize rothermel_surface_fire_spread_no_wind_no_slope
        surface_fire_min        = sf.rothermel_surface_fire_spread_no_wind_no_slope(moisturized_fuel_model)
        residence_time          = surface_fire_min["residence_time"] # min
        reaction_intensity      = surface_fire_min["reaction_intensity"] # Btu/ft^2*min
        # Midflame Wind Speed
        wind_speed_10m          = sqrt(wind_speed_10m_x ** 2.0 + wind_speed_10m_y ** 2.0) # m/min
        wind_speed_20ft         = wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # m/min
        wind_from_direction     = (90.0 - degrees(atan2(wind_speed_10m_y, wind_speed_10m_x)) % 360.0) % 360.0 # deg
        wind_adj_factor         = sf.wind_adjustment_factor(fuel_model["delta"],
                                                            m_to_ft(canopy_height),
                                                            canopy_cover) # unitless
        midflame_wind_speed     = m_to_ft(wind_speed_20ft * wind_adj_factor) # ft/min
        # Max Surface Spread Rate/Direction and Surface Eccentricity
        spread_rate_adjustment  = fuel_spread_adjustment * weather_spread_adjustment # unitless
        surface_fire_max        = sf.rothermel_surface_fire_spread_max(surface_fire_min,
                                                                       midflame_wind_speed,
                                                                       wind_from_direction,
                                                                       slope,
                                                                       aspect,
                                                                       spread_rate_adjustment)
        max_surface_spread_rate = surface_fire_max["max_spread_rate"] # ft/min
        max_spread_direction    = surface_fire_max["max_spread_direction"] # deg
        surface_eccentricity    = surface_fire_max["eccentricity"] # unitless
        # Max Surface Intensity
        #=======================================================================================
        # NOTE: The calculations to determine the fireline_normal_spread_rate have been elided.
        #       Consider ending this function here and making another function to compute the
        #       surface/crown values based on the provided perimeter spread direction.
        #=======================================================================================
        max_flame_depth         = sf.anderson_flame_depth(max_surface_spread_rate, residence_time) # ft
        max_surface_intensity   = Btu_ft_s_to_kW_m(sf.byram_fire_line_intensity(reaction_intensity,
                                                                                max_flame_depth)) # kW/m
        # Check for Crown Fire Initiation
        if cf.van_wagner_crown_fire_initiation(canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture,
                                               max_surface_intensity):
            # Max Crown Spread Rate, Fire Type, and Crown Eccentricity
            max_crown_spread_rate = m_to_ft(cf.cruz_crown_fire_spread(m_min_to_km_hr(wind_speed_10m),
                                                                      canopy_bulk_density,
                                                                      fuel_moisture_dead_1hr)) # ft/min
            fire_type             = 2 if (max_crown_spread_rate < 0.0) else 3 # 2=passive crown, 3=active crown
            max_crown_spread_rate = abs(max_crown_spread_rate) # ft/min
            crown_eccentricity    = (surface_eccentricity
                                     if (max_surface_spread_rate > max_crown_spread_rate)
                                     else crown_fire_eccentricity(m_min_to_mph(wind_speed_20ft))) # unitless
            #=======================================================================================
            # NOTE: The calculations to determine the fireline_normal_spread_rate have been elided.
            #=======================================================================================
            # Max Crown Intensity
            max_crown_intensity   = crown_fire_line_intensity(ft_to_m(max_crown_spread_rate),
                                                              canopy_bulk_density,
                                                              (canopy_height - canopy_base_height),
                                                              Btu_lb_to_kJ_kg(fuel_model["h"][0])) # kW/m
            # Max Combined Spread Rate, Intensity, and Flame Length
            max_spread_rate       = ft_to_m(max(max_surface_spread_rate, max_crown_spread_rate)) # m/min
            max_intensity         = max_surface_intensity + max_crown_intensity # kW/m
            max_flame_length      = ft_to_m(sf.byram_flame_length(kW_m_to_Btu_ft_s(max_intensity))) # m
            # Return Fire Behavior Values
            return {
                "max_spread_rate"        : max_spread_rate,      # m/min
                "max_spread_direction"   : max_spread_direction, # deg
                "max_fire_line_intensity": max_intensity,        # kW/m
                "max_flame_length"       : max_flame_length,     # m
                "fire_type"              : fire_type,            # 2=passive crown, 3=active crown
                "eccentricity"           : crown_eccentricity,   # unitless
            }
        else:
            fire_type                = 1 # 1=surface
            max_surface_flame_length = ft_to_m(sf.byram_flame_length(kW_m_to_Btu_ft_s(max_surface_intensity))) # m
            return {
                "max_spread_rate"        : ft_to_m(max_surface_spread_rate), # m/min
                "max_spread_direction"   : max_spread_direction,             # deg
                "max_fire_line_intensity": max_surface_intensity,            # kW/m
                "max_flame_length"       : max_surface_flame_length,         # m
                "fire_type"              : fire_type,                        # 1=surface
                "eccentricity"           : surface_eccentricity,             # unitless
            }
# burn-cells ends here
