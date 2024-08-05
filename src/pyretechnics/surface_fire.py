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
    - spread_rate (ft/min) orthogonal to the fireline.
    - residence_time (min)
    """
    return spread_rate * residence_time


def calc_fireline_intensity(reaction_intensity, flame_depth):
    """
    Returns the rate of heat release per unit of fire edge in Btu/ft/s given:
    - reaction_intensity (Btu/ft^2/min)
    - flame_depth (ft)
    """
    return (reaction_intensity * flame_depth) / 60.0


def calc_flame_length(fireline_intensity):
    """
    Returns the average flame length in ft given:
    - fireline_intensity (Btu/ft/s)
    """
    return 0.45 * (fireline_intensity ** 0.46)
# surface-fire-intensity-functions ends here
# [[file:../../org/pyretechnics.org::surface-fire-max-effective-wind-speed][surface-fire-max-effective-wind-speed]]
def calc_max_effective_wind_speed(reaction_intensity):
    return 0.9 * reaction_intensity
# surface-fire-max-effective-wind-speed ends here
# [[file:../../org/pyretechnics.org::surface-fire-slope-factor-function][surface-fire-slope-factor-function]]
def get_phi_S_fn(beta):
    if (beta > 0.0):
        G = 5.275 * beta ** -0.3
        return lambda slope: (slope ** 2.0) * G
    else:
        return lambda _: 0.0
# surface-fire-slope-factor-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-wind-factor-function][surface-fire-wind-factor-function]]
def get_phi_W_fn(beta, B, C, F):
    if (beta > 0.0):
        C_over_F = C / F
        return lambda midflame_wind_speed: (midflame_wind_speed ** B) * C_over_F
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
def calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model, spread_rate_adjustment=1.0):
    """
    Given these inputs:
    - moisturized_fuel_model :: dictionary of fuel model and fuel moisture properties
      - delta :: ft                                      (fuel depth)
      - w_o   :: lb/ft^2                                 (ovendry fuel loading)
      - rho_p :: lb/ft^3                                 (ovendry particle density)
      - sigma :: ft^2/ft^3                               (fuel particle surface-area-to-volume ratio)
      - h     :: Btu/lb                                  (fuel particle low heat content)
      - S_T   :: lb minerals/lb ovendry wood             (fuel particle total mineral content)
      - S_e   :: lb silica-free minerals/lb ovendry wood (fuel particle effective mineral content)
      - M_x   :: lb moisture/lb ovendry wood             (fuel particle moisture of extinction)
      - M_f   :: lb moisture/lb ovendry wood             (fuel particle moisture content)
      - f_ij  :: %                                       (percent load per size class)
      - f_i   :: %                                       (percent load per category)
      - g_ij  :: %                                       (percent load per size class - Albini_1976_FIREMOD, page 20)
    - spread_rate_adjustment :: unitless float (1.0 for no adjustment)

    return a dictionary containing these keys:
    - base_spread_rate         :: ft/min
    - base_fireline_intensity  :: Btu/ft/s
    - max_effective_wind_speed :: ft/min
    - get_phi_S                :: lambda: slope => phi_S
    - get_phi_W                :: lambda: midflame_wind_speed => phi_W
    - get_wind_speed           :: lambda: phi_W => midflame_wind_speed
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
    R0             = calc_spread_rate(heat_source, heat_sink)                   # ft/min
    # Calculate base fireline intensity (no wind, no slope)
    t_res          = calc_residence_time(sigma_prime)                           # min
    D_A            = calc_flame_depth(R0, t_res)                                # ft
    I_B            = calc_fireline_intensity(I_R, D_A)                          # Btu/ft/s
    # Pre-compute values related to wind and slope
    U_eff_max      = calc_max_effective_wind_speed(I_R)                         # ft/min
    B              = 0.02526 * (sigma_prime ** 0.54)
    C              = 7.47 * exp(-0.133 * (sigma_prime ** 0.55))
    E              = 0.715 * exp(-3.59 * (sigma_prime / 10000.0))
    F              = (beta / beta_op) ** E
    get_phi_S      = get_phi_S_fn(beta)
    get_phi_W      = get_phi_W_fn(beta, B, C, F)
    get_wind_speed = get_wind_speed_fn(B, C, F)
    # Return no-wind-no-slope surface fire behavior values
    return {
        "base_spread_rate"        : R0 * spread_rate_adjustment,
        "base_fireline_intensity" : I_B * spread_rate_adjustment,
        "max_effective_wind_speed": U_eff_max,
        "get_phi_S"               : get_phi_S,
        "get_phi_W"               : get_phi_W,
        "get_wind_speed"          : get_wind_speed,
    }
# surface-fire-behavior-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::midflame-wind-speed][midflame-wind-speed]]
from math import log, sqrt


def calc_wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover):
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


def calc_midflame_wind_speed(wind_speed_20ft, fuel_bed_depth, canopy_height, canopy_cover):
    """
    wind_speed_20ft :: ft/min
    fuel_bed_depth  :: ft
    canopy_height   :: ft
    canopy_cover    :: 0-1
    """
    wind_adj_factor = calc_wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover)
    return wind_speed_20ft * wind_adj_factor
# midflame-wind-speed ends here
# [[file:../../org/pyretechnics.org::surface-fire-combine-wind-and-slope-vectors][surface-fire-combine-wind-and-slope-vectors]]
from math import radians, degrees, sin, cos, asin, sqrt
from pyretechnics.conversion import opposite_direction


def almost_zero(x):
    return abs(x) < 0.000001


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


def combine_wind_and_slope_vectors(wind_magnitude, downwind_direction, slope_magnitude, upslope_direction):
    difference_angle   = radians((downwind_direction - upslope_direction) % 360.0)
    x_magnitude        = slope_magnitude + wind_magnitude * cos(difference_angle)
    y_magnitude        = wind_magnitude * sin(difference_angle)
    combined_magnitude = sqrt(x_magnitude ** 2.0 + y_magnitude ** 2.0)
    if almost_zero(combined_magnitude):
        return {
            "combined_magnitude": combined_magnitude,
            "combined_direction": upslope_direction,
        }
    else:
        offset             = degrees(asin(abs(y_magnitude) / combined_magnitude))
        offset_prime       = get_offset_prime(x_magnitude, y_magnitude, offset)
        combined_direction = (upslope_direction + offset_prime) % 360.0
        return {
            "combined_magnitude": combined_magnitude,
            "combined_direction": combined_direction,
        }
# surface-fire-combine-wind-and-slope-vectors ends here
# [[file:../../org/pyretechnics.org::surface-fire-eccentricity][surface-fire-eccentricity]]
# FIXME: Surface L/W uses 0.25 but Crown L/W uses 0.125. Check Rothermel 1991.
from math import exp, sqrt
from pyretechnics.conversion import fpm_to_mph


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
# surface-fire-eccentricity ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-max][surface-fire-behavior-max]]
# NOTE: No longer takes ellipse_adjustment_factor parameter
def calc_surface_fire_behavior_max(surface_fire_min, midflame_wind_speed, upwind_direction,
                                   slope, aspect, use_wind_limit=True):
    """
    Given these inputs:
    - surface_fire_min       :: dictionary of no-wind-no-slope surface fire behavior values
      - base_spread_rate         :: ft/min
      - base_fireline_intensity  :: Btu/ft/s
      - max_effective_wind_speed :: ft/min
      - get_phi_S                :: lambda: slope => phi_S
      - get_phi_W                :: lambda: midflame_wind_speed => phi_W
      - get_wind_speed           :: lambda: phi_W => midflame_wind_speed
    - midflame_wind_speed    :: ft/min
    - upwind_direction       :: degrees clockwise from North
    - slope                  :: rise/run
    - aspect                 :: degrees clockwise from North
    - use_wind_limit         :: boolean

    return a dictionary containing these keys:
    - max_spread_rate        :: ft/min
    - max_spread_direction   :: degrees clockwise from North
    - max_fireline_intensity :: Btu/ft/s
    - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
    """
    # Unpack no-wind-no-slope surface fire behavior values
    spread_rate        = surface_fire_min["base_spread_rate"]
    fireline_intensity = surface_fire_min["base_fireline_intensity"]
    max_wind_speed     = surface_fire_min["max_effective_wind_speed"]
    get_phi_W          = surface_fire_min["get_phi_W"]
    get_phi_S          = surface_fire_min["get_phi_S"]
    get_wind_speed     = surface_fire_min["get_wind_speed"]
    # Prepare wind and slope vectors
    phi_W              = get_phi_W(midflame_wind_speed)
    phi_S              = get_phi_S(slope)
    downwind_direction = opposite_direction(upwind_direction)
    upslope_direction  = opposite_direction(aspect)
    # Combine wind and slope vectors
    combined_vector    = combine_wind_and_slope_vectors(phi_W, downwind_direction, phi_S, upslope_direction)
    combined_magnitude = combined_vector["combined_magnitude"]
    combined_direction = combined_vector["combined_direction"]
    # Calculate and return max surface fire behavior values
    effective_wind_speed = get_wind_speed(combined_magnitude)
    if (use_wind_limit and effective_wind_speed > max_wind_speed):
        # Limit effective wind speed to max wind speed
        limited_magnitude    = get_phi_W(max_wind_speed)
        limited_eccentricity = surface_fire_eccentricity(max_wind_speed)
        return {
            "max_spread_rate"       : spread_rate * (1.0 + limited_magnitude),
            "max_spread_direction"  : combined_direction,
            "max_fireline_intensity": fireline_intensity * (1.0 + limited_magnitude),
            "eccentricity"          : limited_eccentricity,
        }
    else:
        combined_eccentricity = surface_fire_eccentricity(effective_wind_speed)
        return {
            "max_spread_rate"       : spread_rate * (1.0 + combined_magnitude),
            "max_spread_direction"  : combined_direction,
            "max_fireline_intensity": fireline_intensity * (1.0 + combined_magnitude),
            "eccentricity"          : combined_eccentricity,
        }
# surface-fire-behavior-max ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-in-direction][surface-fire-behavior-in-direction]]
from math import cos, radians


def smallest_angle_between(theta1, theta2):
  """
  Computes the absolute difference between two angles as an angle between 0° and 180°.
  The return angle has the same cosine as (- theta1 theta2) but may have an opposite sine.
  """
  angle = abs(theta1 - theta2)
  return (360.0 - angle) if (angle > 180.0) else angle


def calc_surface_fire_behavior_in_direction(surface_fire_max, spread_direction):
    # Unpack max surface fire behavior values
    max_spread_rate        = surface_fire_max["max_spread_rate"]
    max_spread_direction   = surface_fire_max["max_spread_direction"]
    max_fireline_intensity = surface_fire_max["max_fireline_intensity"]
    eccentricity           = surface_fire_max["eccentricity"]
    # Calculate adjustment due to the difference angle between spread_direction and max_spread_direction
    theta      = smallest_angle_between(max_spread_direction, spread_direction)
    adjustment = (1.0 - eccentricity) / (1.0 - eccentricity * cos(radians(theta)))
    return {
        "spread_rate"       : max_spread_rate * adjustment,
        "fireline_intensity": max_fireline_intensity * adjustment,
    }
# surface-fire-behavior-in-direction ends here
