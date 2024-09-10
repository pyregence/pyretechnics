# [[file:../../org/pyretechnics.org::van-wagner-critical-fireline-intensity][van-wagner-critical-fireline-intensity]]
def van_wagner_critical_fireline_intensity(canopy_base_height, foliar_moisture):
    """
    Returns the critical fireline intensity (kW/m) given:
    - canopy_base_height :: m
    - foliar_moisture    :: kg moisture/kg ovendry weight

    Constants used:
    460.0 = heat-of-ignition :: kJ/kg
    0.01 = empirical estimate for C in Van Wagner 1977 (eq. 4)
    """
    H = 460.0 + 2600.0 * foliar_moisture
    return (0.01 * canopy_base_height * H) ** 1.5
# van-wagner-critical-fireline-intensity ends here
# [[file:../../org/pyretechnics.org::van-wagner-crown-fire-initiation][van-wagner-crown-fire-initiation]]
def van_wagner_crown_fire_initiation(surface_fireline_intensity, canopy_cover, canopy_base_height, foliar_moisture):
    """
    Returns True if the surface fire transitions to a crown fire or False otherwise given:
    - surface_fireline_intensity :: kW/m
    - canopy_cover               :: 0-1
    - canopy_base_height         :: m
    - foliar_moisture            :: kg moisture/kg ovendry weight
    """
    return (
        surface_fireline_intensity > 0.0
        and
        canopy_cover > 0.4
        and
        surface_fireline_intensity >= van_wagner_critical_fireline_intensity(canopy_base_height, foliar_moisture)
    )
# van-wagner-crown-fire-initiation ends here
# [[file:../../org/pyretechnics.org::cruz-active-crown-fire-spread-rate][cruz-active-crown-fire-spread-rate]]
from math import exp


def cruz_active_crown_fire_spread_rate(wind_speed_10m, canopy_bulk_density, estimated_fine_fuel_moisture):
    """
    Returns the active crown fire spread rate (m/min) given:
    - wind_speed_10m                                   :: km/hr
    - canopy_bulk_density                              :: kg/m^3
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") :: kg moisture/kg ovendry weight
    """
    return (11.02
            * wind_speed_10m ** 0.90
            * canopy_bulk_density ** 0.19
            * exp(-17.0 * estimated_fine_fuel_moisture))
# cruz-active-crown-fire-spread-rate ends here
# [[file:../../org/pyretechnics.org::van-wagner-critical-spread-rate][van-wagner-critical-spread-rate]]
def van_wagner_critical_spread_rate(canopy_bulk_density):
    """
    Returns the critical spread rate (m/min) given:
    - canopy_bulk_density :: kg/m^3
    """
    return 3.0 / canopy_bulk_density
# van-wagner-critical-spread-rate ends here
# [[file:../../org/pyretechnics.org::cruz-passive-crown-fire-spread-rate][cruz-passive-crown-fire-spread-rate]]
from math import exp


def cruz_passive_crown_fire_spread_rate(active_spread_rate, critical_spread_rate):
    """
    Returns the passive crown fire spread rate (m/min) given:
    - active_spread_rate   :: m/min
    - critical_spread_rate :: m/min
    """
    return active_spread_rate * exp(-active_spread_rate / critical_spread_rate)
# cruz-passive-crown-fire-spread-rate ends here
# [[file:../../org/pyretechnics.org::cruz-crown-fire-spread-info][cruz-crown-fire-spread-info]]
def cruz_crown_fire_spread_info(wind_speed_10m, canopy_bulk_density, estimated_fine_fuel_moisture):
    """
    Given these inputs:
    - wind_speed_10m                                   :: km/hr
    - canopy_bulk_density                              :: kg/m^3
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") :: kg moisture/kg ovendry weight

    return a dictionary containing these keys:
    - fire_type            :: "passive_crown" or "active_crown"
    - spread_rate          :: m/min
    - critical_spread_rate :: m/min
    """
    active_spread_rate   = cruz_active_crown_fire_spread_rate(wind_speed_10m,
                                                              canopy_bulk_density,
                                                              estimated_fine_fuel_moisture) # m/min
    critical_spread_rate = van_wagner_critical_spread_rate(canopy_bulk_density) # m/min
    if (active_spread_rate > critical_spread_rate):
        return {
            "fire_type"           : "active_crown",
            "spread_rate"         : active_spread_rate,
            "critical_spread_rate": critical_spread_rate,
        }
    else:
        return {
            "fire_type"           : "passive_crown",
            "spread_rate"         : cruz_passive_crown_fire_spread_rate(active_spread_rate, critical_spread_rate),
            "critical_spread_rate": critical_spread_rate,
        }
# cruz-crown-fire-spread-info ends here
# [[file:../../org/pyretechnics.org::crown-fireline-intensity][crown-fireline-intensity]]
from pyretechnics.surface_fire import calc_flame_length


# NOTE: heat_of_combustion is h from the fuel models (generally 8000 Btu/lb)
# NOTE: ELMFIRE hard-codes heat_of_combustion to 18000 kJ/kg = 7738.6 Btu/lb
def calc_crown_fireline_intensity(crown_spread_rate, canopy_bulk_density, canopy_height,
                                  canopy_base_height, heat_of_combustion):
    """
    Returns the crown fireline intensity (Btu/ft/s OR kW/m) given:
    - crown_spread_rate                                             :: ft/min  OR m/min
    - canopy_bulk_density                                           :: lb/ft^3 OR kg/m^3
    - canopy_height                                                 :: ft      OR m
    - canopy_base_height                                            :: ft      OR m
    - heat_of_combustion                                            :: Btu/lb  OR kJ/kg

    (ft/min * lb/ft^3 * ft * Btu/lb)/60 = (Btu/ft/min)/60 = Btu/ft/s
    OR
    (m/min * kg/m^3 * m * kJ/kg)/60 = (kJ/m*min)/60 = kJ/m*s = kW/m
    """
    canopy_height_difference = canopy_height - canopy_base_height
    return (crown_spread_rate * canopy_bulk_density * canopy_height_difference * heat_of_combustion) / 60.0


def calc_crown_fire_flame_length(surface_fireline_intensity, crown_fireline_intensity):
    """
    Returns the crown fire flame length (m) given:
    - surface_fireline_intensity :: kW/m
    - crown_fireline_intensity   :: kW/m
    """
    return calc_flame_length(surface_fireline_intensity + crown_fireline_intensity) # m
# crown-fireline-intensity ends here
# [[file:../../org/pyretechnics.org::crown-fire-eccentricity][crown-fire-eccentricity]]
from math import sqrt
import pyretechnics.conversion as conv


def crown_length_to_width_ratio(wind_speed_10m, max_length_to_width_ratio=None):
    """
    Calculate the length_to_width_ratio of the crown fire front using eq. 9 from
    Rothermel 1991 given:
    - wind_speed_10m            :: km/hr (aligned with the slope-tangential plane)
    - max_length_to_width_ratio :: float > 0.0 (Optional)
    """
    wind_speed_20ft_mph   = conv.km_hr_to_mph(conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m)) # mph
    length_to_width_ratio = 1.0 + 0.125 * wind_speed_20ft_mph
    if max_length_to_width_ratio:
        return min(length_to_width_ratio, max_length_to_width_ratio)
    else:
        return length_to_width_ratio


def crown_fire_eccentricity(length_to_width_ratio):
    """
    Calculate the eccentricity (E) of the crown fire front using eq. 8 from
    Albini and Chase 1980 given:
    - L/W :: (1: circular spread, > 1: elliptical spread)
    """
    return sqrt(length_to_width_ratio ** 2.0 - 1.0) / length_to_width_ratio
# crown-fire-eccentricity ends here
# [[file:../../org/pyretechnics.org::crown-fire-behavior-max][crown-fire-behavior-max]]
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.surface_fire as sf
import pyretechnics.vector_utils as vu


def calc_crown_fire_behavior_max(canopy_height, canopy_base_height, canopy_bulk_density, heat_of_combustion,
                                 estimated_fine_fuel_moisture, wind_speed_10m, upwind_direction,
                                 slope, aspect, max_length_to_width_ratio=None):
    """
    Given these inputs:
    - canopy_height                                    :: m
    - canopy_base_height                               :: m
    - canopy_bulk_density                              :: kg/m^3
    - heat_of_combustion                               :: kJ/kg
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") :: kg moisture/kg ovendry weight
    - wind_speed_10m                                   :: km/hr
    - upwind_direction                                 :: degrees clockwise from North
    - slope                                            :: rise/run
    - aspect                                           :: degrees clockwise from North
    - max_length_to_width_ratio                        :: float > 0.0 (Optional)

    return a dictionary containing these keys:
    - max_fire_type          :: "passive_crown" or "active_crown"
    - max_spread_rate        :: m/min
    - max_spread_direction   :: (x, y, z) unit vector
    - max_spread_vector      :: (x: m/min, y: m/min, z: m/min)
    - max_fireline_intensity :: kW/m
    - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
    - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
    - critical_spread_rate   :: m/min
    """
    # Reverse the provided wind and slope directions
    downwind_direction = conv.opposite_direction(upwind_direction)
    upslope_direction  = conv.opposite_direction(aspect)
    # Project wind and slope vectors onto the slope-tangential plane
    vectors = sf.project_wind_and_slope_vectors_3d(wind_speed_10m, downwind_direction, slope, upslope_direction)
    wind_vector_3d  = vectors["wind_vector_3d"]  # km/hr
    slope_vector_3d = vectors["slope_vector_3d"] # rise/run
    # Determine the max spread direction
    wind_speed_10m_3d    = vu.vector_magnitude(wind_vector_3d)      # km/hr
    max_spread_direction = (wind_vector_3d / wind_speed_10m_3d      # unit vector in the 3D downwind direction
                            if wind_speed_10m_3d > 0.0
                            else vu.as_unit_vector(slope_vector_3d) # unit vector in the 3D upslope direction
                            if slope > 0.0
                            else np.asarray((0,1,0)))               # default: North
    # Calculate the crown fire behavior in the max spread direction
    spread_info           = cruz_crown_fire_spread_info(wind_speed_10m_3d, canopy_bulk_density,
                                                        estimated_fine_fuel_moisture)
    spread_rate           = spread_info["spread_rate"] # m/min
    fireline_intensity    = calc_crown_fireline_intensity(spread_rate, canopy_bulk_density, canopy_height,
                                                          canopy_base_height, heat_of_combustion) # kW/m
    length_to_width_ratio = crown_length_to_width_ratio(wind_speed_10m_3d, max_length_to_width_ratio) # unitless
    eccentricity          = crown_fire_eccentricity(length_to_width_ratio) # unitless
    return {
        "max_fire_type"         : spread_info["fire_type"],
        "max_spread_rate"       : spread_rate,
        "max_spread_direction"  : max_spread_direction, # unit vector
        "max_spread_vector"     : spread_rate * max_spread_direction,
        "max_fireline_intensity": fireline_intensity,
        "length_to_width_ratio" : length_to_width_ratio,
        "eccentricity"          : eccentricity,
        "critical_spread_rate"  : spread_info["critical_spread_rate"],
    }
# crown-fire-behavior-max ends here
# [[file:../../org/pyretechnics.org::crown-fire-behavior-in-direction][crown-fire-behavior-in-direction]]
import numpy as np


def calc_crown_fire_behavior_in_direction(crown_fire_max, spread_direction):
    """
    Given these inputs:
    - crown_fire_max     :: dictionary of max crown fire behavior values
      - max_fire_type          :: "passive_crown" or "active_crown"
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_spread_vector      :: (x: m/min, y: m/min, z: m/min)
      - max_fireline_intensity :: kW/m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
      - critical_spread_rate   :: m/min
    - spread_direction   :: 3D unit vector on the slope-tangential plane

    return a dictionary containing these keys:
    - fire_type          :: "passive_crown" or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    """
    # Unpack max crown fire behavior values
    max_fire_type          = crown_fire_max["max_fire_type"]
    max_spread_rate        = crown_fire_max["max_spread_rate"]
    max_spread_direction   = crown_fire_max["max_spread_direction"]
    max_fireline_intensity = crown_fire_max["max_fireline_intensity"]
    eccentricity           = crown_fire_max["eccentricity"]
    critical_spread_rate   = crown_fire_max["critical_spread_rate"]
    # Calculate cos(w), where w is the offset angle between these unit vectors on the slope-tangential plane
    cos_w = np.dot(max_spread_direction, np.asarray(spread_direction))
    # Calculate adjustment due to the offset angle from the max spread direction
    adjustment = (1.0 - eccentricity) / (1.0 - eccentricity * cos_w)
    # Adjust the spread rate (possibly switching from an active to passive crown fire)
    spread_rate = max_spread_rate * adjustment
    if spread_rate > critical_spread_rate:
        # Max spread rate was active and directional spread rate remains active
        return {
            "fire_type"         : "active_crown",
            "spread_rate"       : spread_rate,
            "spread_direction"  : spread_direction,
            "fireline_intensity": max_fireline_intensity * adjustment,
        }
    elif max_fire_type == "passive_crown":
        # Max spread rate was passive and directional spread rate remains passive
        return {
            "fire_type"         : "passive_crown",
            "spread_rate"       : spread_rate,
            "spread_direction"  : spread_direction,
            "fireline_intensity": max_fireline_intensity * adjustment,
        }
    else:
        # Max spread rate was active and directional spread rate has become passive
        return {
            "fire_type"         : "passive_crown",
            "spread_rate"       : cruz_passive_crown_fire_spread_rate(spread_rate, critical_spread_rate),
            "spread_direction"  : spread_direction,
            "fireline_intensity": max_fireline_intensity * adjustment,
        }
# crown-fire-behavior-in-direction ends here
# [[file:../../org/pyretechnics.org::combined-fire-behavior][combined-fire-behavior]]
def calc_combined_fire_behavior(surface_fire_behavior, crown_fire_behavior):
    """
    Given these inputs:
    - surface_fire_behavior :: dictionary of surface fire behavior values
      - fire_type              :: "surface"
      - spread_rate            :: m/min
      - spread_direction       :: (x, y, z) unit vector
      - fireline_intensity     :: kW/m
      - flame_length           :: m
    - crown_fire_behavior   :: dictionary of crown fire behavior values
      - fire_type              :: "passive_crown" or "active_crown"
      - spread_rate            :: m/min
      - spread_direction       :: (x, y, z) unit vector
      - fireline_intensity     :: kW/m

    return a dictionary containing these keys:
    - fire_type          :: "surface", "passive_crown", or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    # Unpack the surface fire behavior values
    surface_spread_rate        = surface_fire_behavior["spread_rate"]        # m/min
    surface_spread_direction   = surface_fire_behavior["spread_direction"]   # (x, y, z) unit vector
    surface_fireline_intensity = surface_fire_behavior["fireline_intensity"] # kW/m
    # Unpack the crown fire behavior values
    crown_fire_type          = crown_fire_behavior["fire_type"]          # "passive_crown" or "active_crown"
    crown_spread_rate        = crown_fire_behavior["spread_rate"]        # m/min
    crown_spread_direction   = crown_fire_behavior["spread_direction"]   # (x, y, z) unit vector
    crown_fireline_intensity = crown_fire_behavior["fireline_intensity"] # kW/m
    # Determine whether the surface or crown fire has the fastest spread rate
    if surface_spread_rate > crown_spread_rate:
        # Surface fire spreads faster
        return {
            "fire_type"         : crown_fire_type,
            "spread_rate"       : surface_spread_rate,
            "spread_direction"  : surface_spread_direction,
            "fireline_intensity": surface_fireline_intensity + crown_fireline_intensity,
            "flame_length"      : calc_crown_fire_flame_length(surface_fireline_intensity, crown_fireline_intensity),
        }
    else:
        # Crown fire spreads faster
        return {
            "fire_type"         : crown_fire_type,
            "spread_rate"       : crown_spread_rate,
            "spread_direction"  : crown_spread_direction,
            "fireline_intensity": surface_fireline_intensity + crown_fireline_intensity,
            "flame_length"      : calc_crown_fire_flame_length(surface_fireline_intensity, crown_fireline_intensity),
        }
# combined-fire-behavior ends here
