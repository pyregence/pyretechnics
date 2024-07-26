# [[file:../../org/pyretechnics.org::surface-fire-imports][surface-fire-imports]]
from math import exp
from pyretechnics.fuel_models import map_category, map_size_class, category_sum, size_class_sum
# surface-fire-imports ends here
# [[file:../../org/pyretechnics.org::surface-fire-common-intermediate-calculations][surface-fire-common-intermediate-calculations]]
def calc_surface_area_to_volume_ratio(f_i, f_ij, sigma):
    sigma_prime_i = size_class_sum(lambda i: f_ij[i] * sigma[i])
    return category_sum(lambda i: f_i[i] * sigma_prime_i[i])


def calc_packing_ratio(w_o, rho_p, delta):
    if (delta > 0.0):
        beta_i = size_class_sum(lambda i: w_o[i] / rho_p[i])
        return category_sum(lambda i: beta_i[i]) / delta
    else:
        return 0.0


def calc_optimum_packing_ratio(sigma_prime):
    return (3.348 / sigma_prime ** 0.8189) if (sigma_prime > 0.0) else 1.0
# surface-fire-common-intermediate-calculations ends here
# [[file:../../org/pyretechnics.org::surface-fire-reaction-intensity][surface-fire-reaction-intensity]]
def calc_mineral_damping_coefficients(f_ij, S_e):
    S_e_i = size_class_sum(lambda i: f_ij[i] * S_e[i])
    return map_category(lambda i:
                        (lambda S_e_i:
                         0.174 / (S_e_i ** 0.19) if (S_e_i > 0.0) else 1.0
                         )(S_e_i[i]))


def calc_moisture_damping_coefficients(f_ij, M_f, M_x):
    M_f_i = size_class_sum(lambda i: f_ij[i] * M_f[i])
    M_x_i = size_class_sum(lambda i: f_ij[i] * M_x[i])
    return map_category(lambda i:
                        (lambda M_f, M_x:
                         (lambda r_M:
                          1.0 - 2.59 * r_M + 5.11 * r_M ** 2.0 - 3.52 * r_M ** 3.0
                          )(min(1.0, M_f / M_x) if (M_x > 0.0) else 1.0)
                         )(M_f_i[i], M_x_i[i]))


def calc_low_heat_content(f_ij, h):
    return size_class_sum(lambda i: f_ij[i] * h[i])


def calc_net_fuel_loading(g_ij, w_o, S_T):
    return size_class_sum(lambda i:
                          (lambda g_ij, w_o, S_T:
                           g_ij * w_o * (1.0 - S_T)
                           )(g_ij[i], w_o[i], S_T[i]))


def calc_heat_per_unit_area(eta_S_i, eta_M_i, h_i, W_n_i):
    return category_sum(lambda i: W_n_i[i] * h_i[i] * eta_M_i[i] * eta_S_i[i])


def calc_optimum_reaction_velocity(sigma_prime, beta, beta_op):
    # Albini 1976 replaces 1 / (4.774 * (sigma_prime ** 0.1) - 7.27)
    A               = (133.0 / sigma_prime ** 0.7913) if (sigma_prime > 0.0) else 0.0
    B               = sigma_prime ** 1.5
    C               = beta / beta_op
    # Maximum reaction velocity (1/min)
    Gamma_prime_max = B / (495.0 + 0.0594 * B)
    # Optimum reaction velocity (1/min)
    return Gamma_prime_max * (C ** A) * exp(A * (1.0 - C))


def calc_reaction_intensity(moisturized_fuel_model, sigma_prime, beta, beta_op):
    w_o         = moisturized_fuel_model["w_o"]
    h           = moisturized_fuel_model["h"]
    S_T         = moisturized_fuel_model["S_T"]
    S_e         = moisturized_fuel_model["S_e"]
    M_x         = moisturized_fuel_model["M_x"]
    M_f         = moisturized_fuel_model["M_f"]
    f_ij        = moisturized_fuel_model["f_ij"]
    g_ij        = moisturized_fuel_model["g_ij"]
    eta_S_i     = calc_mineral_damping_coefficients(f_ij, S_e)
    eta_M_i     = calc_moisture_damping_coefficients(f_ij, M_f, M_x)
    h_i         = calc_low_heat_content(f_ij, h)                             # (Btu/lb)
    W_n_i       = calc_net_fuel_loading(g_ij, w_o, S_T)                      # (lb/ft^2)
    Btus        = calc_heat_per_unit_area(eta_S_i, eta_M_i, h_i, W_n_i)      # (Btu/ft^2)
    Gamma_prime = calc_optimum_reaction_velocity(sigma_prime, beta, beta_op) # (1/min)
    return Btus * Gamma_prime                                                # (Btu/ft^2/min)
# surface-fire-reaction-intensity ends here
# [[file:../../org/pyretechnics.org::surface-fire-propagating-flux-ratio][surface-fire-propagating-flux-ratio]]
def calc_propagating_flux_ratio(sigma_prime, beta):
    return exp((0.792 + 0.681 * (sigma_prime ** 0.5)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma_prime)
# surface-fire-propagating-flux-ratio ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-source-no-wind-no-slope][surface-fire-heat-source-no-wind-no-slope]]
def calc_heat_source(I_R, xi):
    return I_R * xi
# surface-fire-heat-source-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::surface-fire-oven-dry-fuel-bed-bulk-density][surface-fire-oven-dry-fuel-bed-bulk-density]]
def calc_ovendry_bulk_density(w_o, delta):
    if (delta > 0.0):
        rho_b_i = size_class_sum(lambda i: w_o[i])
        return category_sum(lambda i: rho_b_i[i]) / delta
    else:
        return 0.0
# surface-fire-oven-dry-fuel-bed-bulk-density ends here
# [[file:../../org/pyretechnics.org::surface-fire-effective-heating-number-distribution][surface-fire-effective-heating-number-distribution]]
def calc_effective_heating_number_distribution(sigma):
    return map_size_class(lambda i:
                          (lambda sigma:
                           exp(-138.0 / sigma) if (sigma > 0.0) else 0.0
                           )(sigma[i]))
# surface-fire-effective-heating-number-distribution ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-of-preignition-distribution][surface-fire-heat-of-preignition-distribution]]
def calc_heat_of_preignition_distribution(M_f):
    return map_size_class(lambda i: 250.0 + 1116.0 * M_f[i])
# surface-fire-heat-of-preignition-distribution ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-sink][surface-fire-heat-sink]]
def calc_heat_sink(f_i, f_ij, rho_b, epsilon_ij, Q_ig_ij):
    effective_heat_of_preignition_i = size_class_sum(lambda i: f_ij[i] * epsilon_ij[i] * Q_ig_ij[i])
    effective_heat_of_preignition   = category_sum(lambda i: f_i[i] * effective_heat_of_preignition_i[i])
    return rho_b * effective_heat_of_preignition
# surface-fire-heat-sink ends here
# [[file:../../org/pyretechnics.org::surface-fire-spread-rate-no-wind-no-slope][surface-fire-spread-rate-no-wind-no-slope]]
def calc_spread_rate(heat_source, heat_sink):
    return heat_source / heat_sink if (heat_sink > 0.0) else 0.0
# surface-fire-spread-rate-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::surface-fire-intensity-functions][surface-fire-intensity-functions]]
def calc_residence_time(sigma_prime):
    return 384.0 / sigma_prime


def calc_flame_depth(spread_rate, residence_time):
    """
    Returns the depth, or front-to-back distance, of the actively flaming zone
    of a free-spreading fire in ft given:
    - spread-rate (ft/min) orthogonal to the fireline.
    - residence-time (min)
    """
    return spread_rate * residence_time


def calc_fire_line_intensity(reaction_intensity, flame_depth):
    """
    Returns the rate of heat release per unit of fire edge in Btu/ft/s given:
    - reaction-intensity (Btu/ft^2/min)
    - flame-depth (ft)
    """
    return (reaction_intensity * flame_depth) / 60.0


def calc_flame_length(fire_line_intensity):
    """
    Returns the average flame length in ft given:
    - fire-line-intensity (Btu/ft/s)
    """
    return 0.45 * (fire_line_intensity ** 0.46)
# surface-fire-intensity-functions ends here
# [[file:../../org/pyretechnics.org::surface-fire-max-effective-wind-speed][surface-fire-max-effective-wind-speed]]
def calc_max_effective_wind_speed(reaction_intensity):
    return 0.9 * reaction_intensity
# surface-fire-max-effective-wind-speed ends here
# [[file:../../org/pyretechnics.org::surface-fire-slope-factor-function][surface-fire-slope-factor-function]]
def get_phi_S_fn(beta):
    if (beta > 0.0):
        G = 5.275 * beta ** -0.3
        return lambda slope: (slope ** 2.0) * G if (slope > 0.0) else 0.0
    else:
        return lambda _: 0.0
# surface-fire-slope-factor-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-wind-factor-function][surface-fire-wind-factor-function]]
def get_phi_W_fn(beta, B, C, F):
    if (beta > 0.0):
        C_over_F = C / F
        return lambda midflame_wind_speed: (midflame_wind_speed ** B) * C_over_F if (midflame_wind_speed > 0.0) else 0.0
    else:
        return lambda _: 0.0
# surface-fire-wind-factor-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-wind-speed-function][surface-fire-wind-speed-function]]
def get_wind_speed_fn(B, C, F):
    F_over_C  = F / C
    B_inverse = 1.0 / B
    return lambda phi_W: (phi_W * F_over_C) ** B_inverse
# surface-fire-wind-speed-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-no-wind-no-slope][surface-fire-behavior-no-wind-no-slope]]
def calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model):
    """
    Given a dictionary containing these keys:
    - delta [fuel depth (ft)]
    - w_o   [ovendry fuel loading (lb/ft^2)]
    - rho_p [ovendry particle density (lb/ft^3)]
    - sigma [fuel particle surface-area-to-volume ratio (ft^2/ft^3)]
    - h     [fuel particle low heat content (Btu/lb)]
    - S_T   [fuel particle total mineral content (lb minerals/lb ovendry wood)]
    - S_e   [fuel particle effective mineral content (lb silica-free minerals/lb ovendry wood)]
    - M_x   [fuel particle moisture of extinction (lb moisture/lb ovendry wood)]
    - M_f   [fuel particle moisture content (lb moisture/lb ovendry wood)]
    - f_ij  [percent of load per size class (%)]
    - f_i   [percent of load per category (%)]
    - g_ij  [percent of load per size class from Albini_1976_FIREMOD, page 20]

    return a dictionary containing these keys:
    - base_spread_rate         (ft/min)
    - base_fireline_intensity  (Btu/ft/s)
    - max_effective_wind_speed (ft/min)
    - get_phi_S                (lambda: slope => phi_S)
    - get_phi_W                (lambda: midflame_wind_speed => phi_W)
    - get_wind_speed           (lambda: phi_W => midflame_wind_speed)
    """
    # Unpack fuel model values
    delta          = moisturized_fuel_model["delta"]
    w_o            = moisturized_fuel_model["w_o"]
    rho_p          = moisturized_fuel_model["rho_p"]
    sigma          = moisturized_fuel_model["sigma"]
    M_f            = moisturized_fuel_model["M_f"]
    f_ij           = moisturized_fuel_model["f_ij"]
    f_i            = moisturized_fuel_model["f_i"]
    # Calculate base spread rate (no wind, no slope)
    sigma_prime    = calc_surface_area_to_volume_ratio(f_i, f_ij, sigma)
    beta           = calc_packing_ratio(w_o, rho_p, delta)
    beta_op        = calc_optimum_packing_ratio(sigma_prime)
    I_R            = calc_reaction_intensity(moisturized_fuel_model, sigma_prime, beta, beta_op) # Btu/ft^2/min
    xi             = calc_propagating_flux_ratio(sigma_prime, beta)
    heat_source    = calc_heat_source(I_R, xi)                                  # Btu/ft^2/min
    rho_b          = calc_ovendry_bulk_density(w_o, delta)                      # lb/ft^3
    epsilon_ij     = calc_effective_heating_number_distribution(sigma)
    Q_ig_ij        = calc_heat_of_preignition_distribution(M_f)                 # Btu/lb
    heat_sink      = calc_heat_sink(f_i, f_ij, rho_b, epsilon_ij, Q_ig_ij)      # Btu/ft^3
    R              = calc_spread_rate(heat_source, heat_sink)                   # ft/min
    # Calculate base fireline intensity (no wind, no slope)
    t_res          = calc_residence_time(sigma_prime)                           # min
    D_A            = calc_flame_depth(R, t_res)                                 # ft
    I_B            = calc_fire_line_intensity(I_R, D_A)                         # Btu/ft/s
    # Pre-compute values related to wind and slope
    U_eff_max      = calc_max_effective_wind_speed(I_R)                         # ft/min
    B              = 0.02526 * (sigma_prime ** 0.54)
    C              = 7.47 * exp(-0.133 * (sigma_prime ** 0.55))
    E              = 0.715 * exp(-3.59 * (sigma_prime / 10000.0))
    F              = (beta / beta_op) ** E
    get_phi_S      = get_phi_S_fn(beta)
    get_phi_W      = get_phi_W_fn(beta, B, C, F)
    get_wind_speed = get_wind_speed_fn(B, C, F)
    return {
        "base_spread_rate"        : R,
        "base_fireline_intensity" : I_B,
        "max_effective_wind_speed": U_eff_max,
        "get_phi_S"               : get_phi_S,
        "get_phi_W"               : get_phi_W,
        "get_wind_speed"          : get_wind_speed,
    }
# surface-fire-behavior-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::wind-adjustment-factor][wind-adjustment-factor]]
from math import log, sqrt

def wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover):
    """
    fuel_bed_depth :: ft
    canopy_height  :: ft
    canopy_cover   :: 0-1
    """
    if (canopy_cover > 0.0) and (canopy_height > 0.0):
        # sheltered: equation 2 based on CC and CH, CR=1 (Andrews 2012)
        A = sqrt((canopy_cover / 3.0) * canopy_height)
        B = log((20.0 + 0.36 * canopy_height) / (0.13 * canopy_height))
        return 0.555 / (A * B)
    elif (fuel_bed_depth > 0.0):
        # unsheltered: equation 6 H_F = H (Andrews 2012)
        A = log((20.0 + 0.36 * fuel_bed_depth) / (0.13 * fuel_bed_depth))
        return 1.83 / A # 1.83 truncated from 1.8328795184533409
    else:
        # non-burnable fuel model
        return 0.0
# wind-adjustment-factor ends here
# [[file:../../org/pyretechnics.org::rothermel-surface-fire-spread-max-and-any][rothermel-surface-fire-spread-max-and-any]]
from math import sin, cos, asin, exp, sqrt, radians, degrees
from pyretechnics.conversion import fpm_to_mph

def almost_zero(x):
    return abs(x) < 0.000001


def smallest_angle_between(theta1, theta2):
  """
  Computes the absolute difference between two angles as an angle between 0° and 180°.
  The return angle has the same cosine as (- theta1 theta2) but may have an opposite sine.
  """
  angle = abs(theta1 - theta2)
  return (360.0 - angle) if (angle > 180.0) else angle


def determine_spread_drivers(midflame_wind_speed, wind_to_direction, slope, slope_direction):
    if almost_zero(slope):
        if almost_zero(midflame_wind_speed):
            return "no_wind_no_slope"
        else:
            return "wind_only"
    elif almost_zero(midflame_wind_speed):
        return "slope_only"
    elif smallest_angle_between(wind_to_direction, slope_direction) < 15.0:
        return "wind_blows_upslope"
    else:
        return "wind_blows_across_slope"


def spread_info_max_no_wind_no_slope(spread_rate):
    return {
        "max_spread_rate"     : spread_rate,
        "max_spread_direction": 0.0,
        "effective_wind_speed": 0.0,
        "eccentricity"        : 0.0,
    }


def spread_info_max_wind_only(spread_rate, phi_W, midflame_wind_speed, wind_to_direction):
    return {
        "max_spread_rate"     : spread_rate * (1.0 + phi_W),
        "max_spread_direction": wind_to_direction,
        "effective_wind_speed": midflame_wind_speed,
        "eccentricity"        : 0.0,
    }


def spread_info_max_slope_only(spread_rate, phi_S, slope_direction, get_wind_speed):
    return {
        "max_spread_rate"     : spread_rate * (1.0 + phi_S),
        "max_spread_direction": slope_direction,
        "effective_wind_speed": get_wind_speed(phi_S),
        "eccentricity"        : 0.0,
    }


def spread_info_max_wind_blows_upslope(spread_rate, phi_combined, slope_direction, get_wind_speed):
    return {
        "max_spread_rate"     : spread_rate * (1.0 + phi_combined),
        "max_spread_direction": slope_direction,
        "effective_wind_speed": get_wind_speed(phi_combined),
        "eccentricity"        : 0.0,
    }


def get_offset_prime(x, y, offset):
    if (x >= 0.0):
        if (y >= 0.0):
            return offset
        else:
            return 360.0 - offset
    elif (y >= 0.0):
        return 180.0 - offset
    else:
        return 180.0 + offset


def spread_info_max_wind_blows_across_slope(spread_rate, phi_W, phi_S, wind_to_direction,
                                            slope_direction, get_wind_speed):
    wind_magnitude     = spread_rate * phi_W
    slope_magnitude    = spread_rate * phi_S
    difference_angle   = radians((wind_to_direction - slope_direction) % 360.0)
    x                  = slope_magnitude + wind_magnitude * cos(difference_angle)
    y                  = wind_magnitude * sin(difference_angle)
    combined_magnitude = sqrt(x * x + y * y)
    if almost_zero(combined_magnitude):
        return {
            "max_spread_rate"     : spread_rate,
            "max_spread_direction": 0.0,
            "effective_wind_speed": 0.0,
            "eccentricity"        : 0.0,
        }
    else:
        max_spread_rate      = spread_rate + combined_magnitude
        phi_combined         = (max_spread_rate / spread_rate) - 1.0
        offset               = degrees(asin(abs(y) / combined_magnitude))
        offset_prime         = get_offset_prime(x, y, offset)
        max_spread_direction = (slope_direction + offset_prime) % 360.0
        effective_wind_speed = get_wind_speed(phi_combined)
        return {
            "max_spread_rate"     : max_spread_rate,
            "max_spread_direction": max_spread_direction,
            "effective_wind_speed": effective_wind_speed,
            "eccentricity"        : 0.0,
        }


def get_spread_info_max(spread_drivers, spread_rate, phi_W, phi_S, midflame_wind_speed,
                        wind_to_direction, slope_direction, get_wind_speed):
    match spread_drivers:
        case "no_wind_no_slope":
            return spread_info_max_no_wind_no_slope(spread_rate)
        case "wind_only":
            return spread_info_max_wind_only(spread_rate, phi_W, midflame_wind_speed, wind_to_direction)
        case "slope_only":
            return spread_info_max_slope_only(spread_rate, phi_S, slope_direction, get_wind_speed)
        case "wind_blows_upslope":
            return spread_info_max_wind_blows_upslope(spread_rate, (phi_W + phi_S), slope_direction, get_wind_speed)
        case "wind_blows_across_slope":
            return spread_info_max_wind_blows_across_slope(spread_rate, phi_W, phi_S, wind_to_direction,
                                                           slope_direction, get_wind_speed)


def scale_spread_to_max_wind_speed(spread_properties, spread_rate, max_wind_speed, phi_max):
    """
    Warning: Mutates spread_properties
    """
    effective_wind_speed = spread_properties["effective_wind_speed"]
    if (effective_wind_speed > max_wind_speed):
        spread_properties["max_spread_rate"]      = spread_rate * (1.0 + phi_max)
        spread_properties["effective_wind_speed"] = max_wind_speed
        return spread_properties
    else:
        return spread_properties


# FIXME: Surface L/W uses 0.25 but Crown L/W uses 0.125. Check Rothermel 1991.
def surface_length_to_width_ratio(effective_wind_speed):
    """
    Calculate the length_to_width_ratio of the surface fire front using eq. 9 from
    Rothermel 1991 given:
    - effective_wind_speed (ft/min)

    L/W = 1 + 0.25 * Ueff_mph
    """
    effective_wind_speed_mph = fpm_to_mph(effective_wind_speed)
    return 1.0 + 0.25 * effective_wind_speed_mph


# FIXME: unused
def surface_length_to_width_ratio_elmfire(effective_wind_speed):
    """
    Calculate the length_to_width_ratio of the surface fire front given:
    - effective_wind_speed (ft/min)

    L/W = min(0.936 * e^(0.2566 * Ueff_mph) + 0.461 * e^(-0.1548 * Ueff_mph) - 0.397, 8.0)
    """
    effective_wind_speed_mph = fpm_to_mph(effective_wind_speed)
    return min((0.936 * exp(0.2566 * effective_wind_speed_mph))
               +
               (0.461 * exp(-0.1548 * effective_wind_speed_mph))
               -
               0.397,
               8.0)


def surface_fire_eccentricity(effective_wind_speed):
    """
    Calculate the eccentricity (E) of the surface fire front using eq. 9 from
    Rothermel 1991 and eq. 8 from Albini and Chase 1980 given:
    - effective_wind_speed (ft/min)

    L/W = 1 + 0.25 * Ueff_mph
    E = sqrt( L/W^2 - 1 ) / L/W
    """
    length_width_ratio = surface_length_to_width_ratio(effective_wind_speed)
    return sqrt(length_width_ratio ** 2.0 - 1.0) / length_width_ratio


def add_eccentricity(spread_properties):
    """
    Warning: Mutates spread_properties
    """
    effective_wind_speed              = spread_properties["effective_wind_speed"] # ft/min
    spread_properties["eccentricity"] = surface_fire_eccentricity(effective_wind_speed)
    return spread_properties


# NOTE: No longer takes ellipse_adjustment_factor parameter
# FIXME: Return max_fireline_intensity
def rothermel_surface_fire_spread_max(surface_fire_min, midflame_wind_speed, wind_from_direction,
                                      slope, aspect, spread_rate_adjustment=1.0):
    spread_rate        = surface_fire_min["base_spread_rate"] * spread_rate_adjustment
    fireline_intensity = surface_fire_min["base_fireline_intensity"]
    max_wind_speed     = surface_fire_min["max_effective_wind_speed"]
    get_phi_S          = surface_fire_min["get_phi_S"]
    get_phi_W          = surface_fire_min["get_phi_W"]
    get_wind_speed     = surface_fire_min["get_wind_speed"]
    slope_direction    = (aspect + 180.0) % 360.0
    wind_to_direction  = (wind_from_direction + 180.0) % 360.0
    phi_S              = get_phi_S(slope)
    phi_W              = get_phi_W(midflame_wind_speed)
    phi_max            = get_phi_W(max_wind_speed)
    spread_drivers     = determine_spread_drivers(midflame_wind_speed, wind_to_direction, slope, slope_direction)
    spread_info_max    = get_spread_info_max(spread_drivers, spread_rate, phi_W, phi_S, midflame_wind_speed,
                                             wind_to_direction, slope_direction, get_wind_speed)
    return add_eccentricity(scale_spread_to_max_wind_speed(spread_info_max, spread_rate, max_wind_speed, phi_max))


def compute_spread_rate(max_spread_rate, max_spread_direction, eccentricity, spread_direction):
    theta = smallest_angle_between(max_spread_direction, spread_direction)
    if almost_zero(eccentricity) or almost_zero(theta):
        return max_spread_rate
    else:
        return max_spread_rate * (1.0 - eccentricity) / (1.0 - eccentricity * cos(radians(theta)))
# rothermel-surface-fire-spread-max-and-any ends here
