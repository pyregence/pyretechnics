# [[file:../../org/pyretechnics.org::van-wagner-crown-fire-initiation][van-wagner-crown-fire-initiation]]
def van_wagner_critical_fire_line_intensity(canopy_base_height, foliar_moisture):
    """
    Ouputs the critical fire line intensity (kW/m) using:
    - canopy_base_height (m)
    - foliar_moisture (0-1)

    Constants used:
    460.0 = heat-of-ignition (kJ/kg)
    0.01 = empirical estimate for C in eq. 4
    """
    return ((foliar_moisture * 2600.0 + 460.0) * 0.01 * canopy_base_height) ** 1.5


def van_wagner_crown_fire_initiation(canopy_cover, canopy_base_height, foliar_moisture, fire_line_intensity):
    """
    - canopy_cover (0-1)
    - canopy_base_height (m)
    - foliar_moisture (0-1)
    - fire_line_intensity (kW/m)
    """
    return ((canopy_cover > 0.4)
            and
            (fire_line_intensity > 0.0)
            and
            (canopy_base_height > 0.0)
            and
            (fire_line_intensity >= van_wagner_critical_fire_line_intensity(canopy_base_height, foliar_moisture)))
# van-wagner-crown-fire-initiation ends here
# [[file:../../org/pyretechnics.org::cruz-crown-fire-spread][cruz-crown-fire-spread]]
from math import exp, prod

def cruz_active_crown_fire_spread(wind_speed_10m, canopy_bulk_density, estimated_fine_fuel_moisture):
    """
    Returns active spread-rate in m/min given:
    - wind_speed_10m (km/hr)
    - canopy_bulk_density (kg/m^3)
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") (0-1)
    """
    return prod(11.02,
                wind_speed_10m ** 0.90,
                canopy_bulk_density ** 0.19,
                exp(-17.0 * estimated_fine_fuel_moisture))


def cruz_passive_crown_fire_spread(active_spread_rate, critical_spread_rate):
    """
    Returns passive spread-rate in m/min given:
    - active_spread_rate (m/min)
    - critical_spread_rate (m/min)
    """
    return active_spread_rate * exp(-1.0 * active_spread_rate / critical_spread_rate)


def cruz_crown_fire_spread(wind_speed_10m, canopy_bulk_density, estimated_fine_fuel_moisture):
    """
    Returns spread-rate in m/min given:
    - wind_speed_10m (km/hr)
    - canopy_bulk_density (kg/m^3)
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") (0-1)
    NOTE: A positive spread-rate indicates active crowning.
          A negative spread-rate indicates passive crowning.
    """
    active_spread_rate   = cruz_active_crown_fire_spread(wind_speed_10m,
                                                         canopy_bulk_density,
                                                         estimated_fine_fuel_moisture)
    critical_spread_rate = 3.0 / canopy_bulk_density # m/min
    if (active_spread_rate > critical_spread_rate):
        return active_spread_rate
    else:
        # NOTE: Use minus as passive flag
        return -1.0 * cruz_passive_crown_fire_spread(active_spread_rate, critical_spread_rate)
# cruz-crown-fire-spread ends here
# [[file:../../org/pyretechnics.org::crown-fire-line-intensity][crown-fire-line-intensity]]
# NOTE: heat_of_combustion is h from the fuel models (generally 8000 Btu/lb)
# NOTE: ELMFIRE hard-codes heat_of_combustion to 18000 kJ/kg = 7738.6 Btu/lb
def crown_fire_line_intensity(crown_spread_rate, canopy_bulk_density, canopy_height_difference, heat_of_combustion):
    """
    Returns the crown_fire_line_intensity in Btu/ft*s OR kW/m, given:
    - crown_spread_rate (ft/min OR m/min)
    - canopy_bulk_density (lb/ft^3 OR kg/m^3)
    - canopy_height_difference (canopy_height - canopy_base_height) (ft OR m)
    - heat_of_combustion (Btu/lb OR kJ/kg)

    (ft/min * lb/ft^3 * ft * Btu/lb)/60 = (Btu/ft*min)/60 = Btu/ft*s
    OR
    (m/min * kg/m^3 * m * kJ/kg)/60 = (kJ/m*min)/60 = kJ/m*s = kW/m
    """
    return (crown_spread_rate * canopy_bulk_density * canopy_height_difference * heat_of_combustion) / 60.0
# crown-fire-line-intensity ends here
# [[file:../../org/pyretechnics.org::crown-eccentricity][crown-eccentricity]]
from math import sqrt

# NOTE: No longer takes ellipse_adjustment_factor argument
# FIXME: Surface L/W uses 0.25 but Crown L/W uses 0.125. Check Rothermel 1991.
def crown_length_to_width_ratio(wind_speed_20ft):
    """
    Calculate the length_to_width_ratio of the crown fire front using eq. 9 from
    Rothermel 1991 given:
    - wind_speed_20ft (mph)

    L/W = 1 + 0.125 * U20_mph
    """
    return 1.0 + 0.125 * wind_speed_20ft


# FIXME: unused
def crown_length_to_width_ratio_elmfire(wind_speed_20ft, max_length_to_width_ratio):
    """
    Calculate the length_to_width_ratio of the crown fire front using eq. 9 from
    Rothermel 1991 given:
    - wind_speed_20ft (mph)
    - max_length_to_width_ratio (int > 0)

    L/W = min(1.0 + 0.125 * U20_mph, L/W_max)
    """
    return min((1.0 + 0.125 * wind_speed_20ft), max_length_to_width_ratio)


# NOTE: No longer takes ellipse_adjustment_factor argument
def crown_fire_eccentricity(wind_speed_20ft):
    """
    Calculate the eccentricity (E) of the crown fire front using eq. 9 from
    Rothermel 1991 and eq. 8 from Albini and Chase 1980 given:
    - wind_speed_20ft (mph)

    L/W = 1 + 0.125 * U20_mph
    E = sqrt( L/W^2 - 1 ) / L/W
    """
    length_width_ratio = crown_length_to_width_ratio(wind_speed_20ft)
    return sqrt(length_width_ratio ** 2.0 - 1.0) / length_width_ratio
# crown-eccentricity ends here
