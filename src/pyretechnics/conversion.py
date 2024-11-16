# [[file:../../org/pyretechnics.org::units-conversion][units-conversion]]
# TODO: Fix error with importing pyretechnics.types
import cython
if cython.compiled:
    from cython.cimports.pyretechnics.math import sqrt, sin, cos, tan, atan, atan2
    from cython.cimports.pyretechnics.types import vec_xy
else:
    from math import sqrt, sin, cos, tan, atan, atan2
    from pyretechnics.types import vec_xy


from math import degrees, radians, pi
import cython as cy


def F_to_K(degrees):
    """Convert fahrenheit to kelvin."""
    return (degrees + 459.67) * 0.5555555555555556


def K_to_F(degrees):
    """Convert kelvin to fahrenheit."""
    return (degrees * 1.8) - 459.67


def F_to_C(degrees):
    """Convert fahrenheit to celsius."""
    return (degrees - 32.0) * 0.5555555555555556


def C_to_F(degrees):
    """Convert celsius to fahrenheit."""
    return (degrees * 1.8) + 32.0


def deg_to_ratio(degrees):
    """Convert degrees to ratio."""
    return tan(radians(degrees))


def ratio_to_deg(ratio):
    """Convert ratio to degrees."""
    return degrees(atan(ratio))


def ch_to_m(ch):
    """Convert chains to meters."""
    return ch * 20.1168


def m_to_ch(m):
    """Convert meters to chains."""
    return m * 0.0497097


def m_to_ft(m):
    """Convert meters to feet."""
    return m * 3.281


def ft_to_m(ft):
    """Convert feet to meters."""
    return ft * 0.30478512648582745


def mph_to_mps(mph):
    """Convert miles per hour to meters per second."""
    return mph * 0.44701818551254696


def mps_to_mph(mps):
    """Convert meters per second to miles per hour."""
    return mps * 2.237045454545455


def km_hr_to_mps(km_hr):
    """Convert kilometers per hour to meters per second."""
    return km_hr * 0.277764222883701


def mps_to_km_hr(mps):
    """Convert meters per second to kilometers per hour."""
    return mps * 3.6001756800000004


def mph_to_km_hr(mph):
    """Convert miles per hour to kilometers per hour."""
    return mph * 1.609344


def km_hr_to_mph(km_hr):
    """Convert kilometers per hour to miles per hour."""
    return km_hr * 0.621371192237334


def m_min_to_km_hr(m_min):
    """Convert meters per minute to kilometers per hour."""
    return m_min * 0.06


def km_hr_to_m_min(km_hr):
    """Convert kilometers per hour to meters per minute."""
    return km_hr / 0.06


def m_min_to_mph(m_min):
    """Convert meters per minute to miles per hour."""
    return m_min * 0.0372840909091


def mph_to_m_min(mph):
    """Convert miles per hour to meters per minute."""
    return mph * 26.8210911307


def mps_to_fpm(mps):
    """Convert meters per second to feet per minute."""
    return mps * 196.86


def fpm_to_mps(fpm):
    """Convert feet per minute to meters per second."""
    return fpm / 196.86


def mph_to_fpm(mph):
    """Convert miles per hour to feet per minute."""
    return mph * 88.0


def fpm_to_mph(fpm):
    """Convert feet per minute to miles per hour."""
    return fpm / 88.0


def Btu_ft_s_to_kW_m(Btu_ft_s):
    """Convert BTU per feet per second to kilowatt per meter."""
    return Btu_ft_s * 3.46165186


def kW_m_to_Btu_ft_s(kW_m):
    """Convert kilowatt per meter to BTU per feet per second."""
    return kW_m * 0.28887942532730604


def Btu_lb_to_kJ_kg(Btu_lb):
    """Convert BTU per lb to kilojoule per kilogram."""
    return Btu_lb * 2.3259999996185


def kJ_kg_to_Btu_lb(kJ_kg):
    """Convert kilojoule per kilogram to BTU per lb."""
    return kJ_kg / 2.3259999996185


def kg_m3_to_lb_ft3(kg_m3):
    """Convert kilogram per cubic meter to pound per cubic foot."""
    return kg_m3 * 0.0624


def lb_ft3_to_kg_m3(lb_ft3):
    """Convert pound per cubic foot to kilogram per cubic meter."""
    return lb_ft3 * 16.025641025641026


def percent_to_dec(percent):
    """Convert percent to decimal."""
    return percent * 0.01


def dec_to_percent(decimal):
    """Convert decimal to percent."""
    return decimal * 100.0


def sec_to_min(seconds):
    """Convert seconds to minutes."""
    return seconds * 0.016666666666666666


def min_to_sec(minutes):
    """Convert minutes to seconds."""
    return minutes * 60.0


def ms_to_min(milliseconds):
    """Convert milliseconds to minutes."""
    return milliseconds * 0.000016667


def min_to_ms(minutes):
    """Convert minutes to milliseconds."""
    return int(minutes * 60000.0)


def hour_to_min(hours):
    """Converts hours to minutes."""
    return hours * 60.0


def min_to_hour(minutes):
    """Converts minutes to hours. (rounds down)"""
    return int(minutes / 60.0)


def day_to_min(days):
    """Convert days to minutes."""
    return days * 1440.0


def min_to_day(minutes):
    """Convert minutes to days."""
    return minutes / 1440.0


def cartesian_to_polar(x, y):
    """Convert cartesian coordinates (x, y) to polar coordinates (r, theta)."""
    r         = sqrt(x ** 2.0 + y ** 2.0)
    theta_rad = atan2(y, x)
    theta     = degrees(theta_rad) % 360.0
    return (r, theta)


def polar_to_cartesian(r, theta):
    """Convert polar coordinates (r, theta) to cartesian coordinates (x, y)."""
    theta_rad = radians(theta)
    x         = r * cos(theta_rad)
    y         = r * sin(theta_rad)
    return (x, y)


# TODO: Implement a primitive degrees function
@cy.profile(True)
@cy.ccall
def cartesian_to_azimuthal(x: cy.float, y: cy.float) -> vec_xy:
    """Convert cartesian coordinates (x, y) to azimuthal coordinates (r, azimuth)."""
    r          : cy.float = sqrt(x * x + y * y)
    azimuth_rad: cy.float = atan2(x, y)
    azimuth    : cy.float = degrees(azimuth_rad) % 360.0
    return (r, azimuth)


# TODO: Implement a primitive radians function
@cy.profile(True)
@cy.ccall
def azimuthal_to_cartesian(r: cy.float, azimuth: cy.float) -> vec_xy:
    """Convert azimuthal coordinates (r, azimuth) to cartesian coordinates (x, y)."""
    azimuth_rad: cy.float = radians(azimuth)
    x          : cy.float = r * sin(azimuth_rad)
    y          : cy.float = r * cos(azimuth_rad)
    return (x, y)


@cy.profile(True)
@cy.ccall
def opposite_direction(theta: cy.float) -> cy.float:
    """Convert theta to theta + 180 degrees."""
    return (theta + 180.0) % 360.0


def wind_speed_10m_to_wind_speed_20ft(wind_speed_10m):
    """Convert wind speed at 10m to wind speed at 20ft."""
    return 0.87 * wind_speed_10m


def wind_speed_20ft_to_wind_speed_10m(wind_speed_20ft):
    """Convert wind speed at 20ft to wind speed at 10m."""
    return wind_speed_20ft / 0.87
# units-conversion ends here
