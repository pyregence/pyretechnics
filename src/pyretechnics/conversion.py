# [[file:../../org/pyretechnics.org::units-conversion][units-conversion]]
import cython
import cython as cy
if cython.compiled:
    from cython.cimports.libc.math import pi, sqrt, sin, cos, tan, atan, atan2
    from cython.cimports.pyretechnics.cy_types import vec_xy
else:
    from math import pi, sqrt, sin, cos, tan, atan, atan2
    from pyretechnics.py_types import vec_xy


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def rad_to_deg(radians: cy.float) -> cy.float:
    """Convert radians to degrees."""
    return radians * 180.0 / pi


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def deg_to_rad(degrees: cy.float) -> cy.float:
    """Convert degrees to radians."""
    return degrees * pi / 180.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def deg_to_ratio(degrees: cy.float) -> cy.float:
    """Convert degrees to ratio."""
    return tan(deg_to_rad(degrees))


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def ratio_to_deg(ratio: cy.float) -> cy.float:
    """Convert ratio to degrees."""
    return rad_to_deg(atan(ratio))


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def F_to_K(degrees: cy.float) -> cy.float:
    """Convert fahrenheit to kelvin."""
    return (degrees + 459.67) * 0.5555555555555556


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def K_to_F(degrees: cy.float) -> cy.float:
    """Convert kelvin to fahrenheit."""
    return (degrees * 1.8) - 459.67


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def F_to_C(degrees: cy.float) -> cy.float:
    """Convert fahrenheit to celsius."""
    return (degrees - 32.0) * 0.5555555555555556


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def C_to_F(degrees: cy.float) -> cy.float:
    """Convert celsius to fahrenheit."""
    return (degrees * 1.8) + 32.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def ch_to_m(ch: cy.float) -> cy.float:
    """Convert chains to meters."""
    return ch * 20.1168


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def m_to_ch(m: cy.float) -> cy.float:
    """Convert meters to chains."""
    return m * 0.0497097


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def m_to_ft(m: cy.float) -> cy.float:
    """Convert meters to feet."""
    return m * 3.281


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def ft_to_m(ft: cy.float) -> cy.float:
    """Convert feet to meters."""
    return ft * 0.30478512648582745


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mph_to_mps(mph: cy.float) -> cy.float:
    """Convert miles per hour to meters per second."""
    return mph * 0.44701818551254696


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mps_to_mph(mps: cy.float) -> cy.float:
    """Convert meters per second to miles per hour."""
    return mps * 2.237045454545455


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def km_hr_to_mps(km_hr: cy.float) -> cy.float:
    """Convert kilometers per hour to meters per second."""
    return km_hr * 0.277764222883701


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mps_to_km_hr(mps: cy.float) -> cy.float:
    """Convert meters per second to kilometers per hour."""
    return mps * 3.6001756800000004


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mph_to_km_hr(mph: cy.float) -> cy.float:
    """Convert miles per hour to kilometers per hour."""
    return mph * 1.609344


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def km_hr_to_mph(km_hr: cy.float) -> cy.float:
    """Convert kilometers per hour to miles per hour."""
    return km_hr * 0.621371192237334


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def m_min_to_km_hr(m_min: cy.float) -> cy.float:
    """Convert meters per minute to kilometers per hour."""
    return m_min * 0.06


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def km_hr_to_m_min(km_hr: cy.float) -> cy.float:
    """Convert kilometers per hour to meters per minute."""
    return km_hr / 0.06


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def m_min_to_mph(m_min: cy.float) -> cy.float:
    """Convert meters per minute to miles per hour."""
    return m_min * 0.0372840909091


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mph_to_m_min(mph: cy.float) -> cy.float:
    """Convert miles per hour to meters per minute."""
    return mph * 26.8210911307


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mps_to_fpm(mps: cy.float) -> cy.float:
    """Convert meters per second to feet per minute."""
    return mps * 196.86


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def fpm_to_mps(fpm: cy.float) -> cy.float:
    """Convert feet per minute to meters per second."""
    return fpm / 196.86


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def mph_to_fpm(mph: cy.float) -> cy.float:
    """Convert miles per hour to feet per minute."""
    return mph * 88.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def fpm_to_mph(fpm: cy.float) -> cy.float:
    """Convert feet per minute to miles per hour."""
    return fpm / 88.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def Btu_ft_s_to_kW_m(Btu_ft_s: cy.float) -> cy.float:
    """Convert BTU per feet per second to kilowatt per meter."""
    return Btu_ft_s * 3.46165186


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def kW_m_to_Btu_ft_s(kW_m: cy.float) -> cy.float:
    """Convert kilowatt per meter to BTU per feet per second."""
    return kW_m * 0.28887942532730604


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def Btu_lb_to_kJ_kg(Btu_lb: cy.float) -> cy.float:
    """Convert BTU per lb to kilojoule per kilogram."""
    return Btu_lb * 2.3259999996185


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def kJ_kg_to_Btu_lb(kJ_kg: cy.float) -> cy.float:
    """Convert kilojoule per kilogram to BTU per lb."""
    return kJ_kg / 2.3259999996185


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def kg_m3_to_lb_ft3(kg_m3: cy.float) -> cy.float:
    """Convert kilogram per cubic meter to pound per cubic foot."""
    return kg_m3 * 0.0624


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def lb_ft3_to_kg_m3(lb_ft3: cy.float) -> cy.float:
    """Convert pound per cubic foot to kilogram per cubic meter."""
    return lb_ft3 * 16.025641025641026


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def percent_to_dec(percent: cy.float) -> cy.float:
    """Convert percent to decimal."""
    return percent * 0.01


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def dec_to_percent(decimal: cy.float) -> cy.float:
    """Convert decimal to percent."""
    return decimal * 100.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def sec_to_min(seconds: cy.float) -> cy.float:
    """Convert seconds to minutes."""
    return seconds * 0.016666666666666666


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def min_to_sec(minutes: cy.float) -> cy.float:
    """Convert minutes to seconds."""
    return minutes * 60.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def ms_to_min(milliseconds: cy.float) -> cy.float:
    """Convert milliseconds to minutes."""
    return milliseconds * 0.000016667


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def min_to_ms(minutes: cy.float) -> cy.float:
    """Convert minutes to milliseconds."""
    return minutes * 60000.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def hour_to_min(hours: cy.float) -> cy.float:
    """Converts hours to minutes."""
    return hours * 60.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def min_to_hour(minutes: cy.float) -> cy.float:
    """Converts minutes to hours."""
    return minutes / 60.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def day_to_min(days: cy.float) -> cy.float:
    """Convert days to minutes."""
    return days * 1440.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def min_to_day(minutes: cy.float) -> cy.float:
    """Convert minutes to days."""
    return minutes / 1440.0


# TODO Return vec_rt
@cy.ccall
@cy.exceptval(check=False)
def cartesian_to_polar(x: cy.float, y: cy.float) -> vec_xy:
    """Convert cartesian coordinates (x, y) to polar coordinates (r, theta)."""
    r        : cy.float = sqrt(x * x + y * y)
    theta_rad: cy.float = atan2(y, x)
    theta    : cy.float = (rad_to_deg(theta_rad) + 360.0) % 360.0
    return (r, theta)


@cy.ccall
@cy.exceptval(check=False)
def polar_to_cartesian(r: cy.float, theta: cy.float) -> vec_xy:
    """Convert polar coordinates (r, theta) to cartesian coordinates (x, y)."""
    theta_rad: cy.float = deg_to_rad(theta)
    x        : cy.float = r * cos(theta_rad)
    y        : cy.float = r * sin(theta_rad)
    return (x, y)


# TODO Return vec_ra
@cy.ccall
@cy.exceptval(check=False)
def cartesian_to_azimuthal(x: cy.float, y: cy.float) -> vec_xy:
    """Convert cartesian coordinates (x, y) to azimuthal coordinates (r, azimuth)."""
    r          : cy.float = sqrt(x * x + y * y)
    azimuth_rad: cy.float = atan2(x, y)
    azimuth    : cy.float = (rad_to_deg(azimuth_rad) + 360.0) % 360.0
    return (r, azimuth)


@cy.ccall
@cy.exceptval(check=False)
def azimuthal_to_cartesian(r: cy.float, azimuth: cy.float) -> vec_xy:
    """Convert azimuthal coordinates (r, azimuth) to cartesian coordinates (x, y)."""
    azimuth_rad: cy.float = deg_to_rad(azimuth)
    x          : cy.float = r * sin(azimuth_rad)
    y          : cy.float = r * cos(azimuth_rad)
    return (x, y)


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def opposite_direction(theta: cy.float) -> cy.float:
    """Convert theta to theta + 180 degrees."""
    return (theta + 180.0) % 360.0


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def wind_speed_10m_to_wind_speed_20ft(wind_speed_10m: cy.float) -> cy.float:
    """Convert wind speed at 10m to wind speed at 20ft."""
    return wind_speed_10m / 1.15


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def wind_speed_20ft_to_wind_speed_10m(wind_speed_20ft: cy.float) -> cy.float:
    """Convert wind speed at 20ft to wind speed at 10m."""
    return wind_speed_20ft * 1.15
# units-conversion ends here
