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
import pyretechnics.conversion as conv


def calc_residence_time(sigma_prime):
    return 384.0 / sigma_prime if (sigma_prime > 0.0) else 0.0


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
    Returns the average flame length (m) given:
    - fireline_intensity :: kW/m
    """
    return conv.ft_to_m(0.45 * (conv.kW_m_to_Btu_ft_s(fireline_intensity) ** 0.46))
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
import pyretechnics.conversion as conv


def get_phi_W_fn(B, C, F):
    if (F > 0.0):
        C_over_F = C / F
        return lambda midflame_wind_speed: (conv.m_to_ft(midflame_wind_speed) ** B) * C_over_F
    else:
        return lambda _: 0.0
# surface-fire-wind-factor-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-wind-speed-function][surface-fire-wind-speed-function]]
import pyretechnics.conversion as conv


def get_wind_speed_fn(B, C, F):
    if (B > 0.0):
        F_over_C  = F / C
        B_inverse = 1.0 / B
        return lambda phi_W: conv.ft_to_m((phi_W * F_over_C) ** B_inverse)
    else:
        return lambda _: 0.0
# surface-fire-wind-speed-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-no-wind-no-slope][surface-fire-behavior-no-wind-no-slope]]
import pyretechnics.conversion as conv


def calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model, spread_rate_adjustment=1.0):
    """
    Given these inputs:
    - moisturized_fuel_model :: dictionary of fuel model and fuel moisture properties
      - delta :: ft                                        (fuel depth)
      - w_o   :: lb/ft^2                                   (ovendry fuel loading)
      - rho_p :: lb/ft^3                                   (ovendry particle density)
      - sigma :: ft^2/ft^3                                 (fuel particle surface-area-to-volume ratio)
      - h     :: Btu/lb                                    (fuel particle low heat content)
      - S_T   :: lb minerals/lb ovendry weight             (fuel particle total mineral content)
      - S_e   :: lb silica-free minerals/lb ovendry weight (fuel particle effective mineral content)
      - M_x   :: lb moisture/lb ovendry weight             (fuel particle moisture of extinction)
      - M_f   :: lb moisture/lb ovendry weight             (fuel particle moisture content)
      - f_ij  :: %                                         (percent load per size class)
      - f_i   :: %                                         (percent load per category)
      - g_ij  :: %                                         (percent load per size class - Albini_1976_FIREMOD, page 20)
    - spread_rate_adjustment :: unitless float (1.0 for no adjustment)

    return a dictionary containing these keys:
    - base_spread_rate         :: m/min
    - base_fireline_intensity  :: kW/m
    - max_effective_wind_speed :: m/min
    - get_phi_S                :: lambda: slope (rise/run) => phi_S (unitless)
    - get_phi_W                :: lambda: midflame_wind_speed (m/min) => phi_W (unitless)
    - get_wind_speed           :: lambda: phi_W (unitless) => midflame_wind_speed (m/min)
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
    I_s            = calc_fireline_intensity(I_R, D_A)                          # Btu/ft/s
    # Pre-compute values related to wind and slope
    U_eff_max      = calc_max_effective_wind_speed(I_R)                         # ft/min
    B              = 0.02526 * (sigma_prime ** 0.54)
    C              = 7.47 * exp(-0.133 * (sigma_prime ** 0.55))
    E              = 0.715 * exp(-3.59 * (sigma_prime / 10000.0))
    F              = (beta / beta_op) ** E
    get_phi_S      = get_phi_S_fn(beta)
    get_phi_W      = get_phi_W_fn(B, C, F)
    get_wind_speed = get_wind_speed_fn(B, C, F)
    # Return no-wind-no-slope surface fire behavior values
    return {
        "base_spread_rate"        : conv.ft_to_m(R0 * spread_rate_adjustment),
        "base_fireline_intensity" : conv.Btu_ft_s_to_kW_m(I_s * spread_rate_adjustment),
        "max_effective_wind_speed": conv.ft_to_m(U_eff_max),
        "get_phi_S"               : get_phi_S,
        "get_phi_W"               : get_phi_W,
        "get_wind_speed"          : get_wind_speed,
    }
# surface-fire-behavior-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::midflame-wind-speed][midflame-wind-speed]]
from math import log, sqrt


def calc_wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover):
    """
    Return the wind adjustment factor (unitless) given these inputs:
    - fuel_bed_depth :: ft
    - canopy_height  :: ft
    - canopy_cover   :: 0-1
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
    Return the midflame wind speed (S) given these inputs:
    - wind_speed_20ft :: S
    - fuel_bed_depth  :: ft
    - canopy_height   :: ft
    - canopy_cover    :: 0-1
    """
    wind_adj_factor = calc_wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover)
    return wind_speed_20ft * wind_adj_factor
# midflame-wind-speed ends here
# [[file:../../org/pyretechnics.org::surface-fire-combine-wind-and-slope-vectors][surface-fire-combine-wind-and-slope-vectors]]
import numpy as np
from pyretechnics.conversion import azimuthal_to_cartesian
from pyretechnics.vector_utils import vector_magnitude, as_unit_vector, to_slope_plane


def project_wind_and_slope_vectors_3d(wind_speed, downwind_direction, slope, upslope_direction):
    """
    Given these inputs:
    - wind_speed         :: S
    - downwind_direction :: degrees clockwise from North
    - slope              :: rise/run
    - upslope_direction  :: degrees clockwise from North

    return a dictionary containing these keys:
    - wind_vector_3d  :: (x: S, y: S, z: S)
    - slope_vector_3d :: (x, y, z)
    """
    # Convert wind and slope vectors from azimuthal to cartesian coordinates
    wind_vector_2d  = azimuthal_to_cartesian(wind_speed, downwind_direction)
    slope_vector_2d = azimuthal_to_cartesian(slope, upslope_direction)
    # Project wind and slope vectors onto the slope-tangential plane
    wind_vector_3d  = to_slope_plane(wind_vector_2d, slope_vector_2d)
    slope_vector_3d = to_slope_plane(slope_vector_2d, slope_vector_2d)
    return {
        "wind_vector_3d" : wind_vector_3d,
        "slope_vector_3d": slope_vector_3d,
    }


def get_phi_E(wind_vector_3d, slope_vector_3d, phi_W, phi_S):
    # Convert wind and slope vectors to unit vectors on the slope-tangential plane
    w_S = as_unit_vector(wind_vector_3d)  if phi_W > 0.0 else wind_vector_3d
    u_S = as_unit_vector(slope_vector_3d) if phi_S > 0.0 else slope_vector_3d
    # Create the 3D slope-tangential phi_W, phi_S, and phi_E vectors
    phi_W_3d = phi_W * w_S
    phi_S_3d = phi_S * u_S
    phi_E_3d = phi_W_3d + phi_S_3d
    # Calculate phi_E
    phi_E = vector_magnitude(phi_E_3d)
    # Determine max spread direction and return results
    if phi_E > 0.0:
        return {
            "phi_E"               : phi_E,
            "max_spread_direction": phi_E_3d / phi_E,
        }
    elif phi_S > 0.0:
        return {
            "phi_E"               : phi_E,
            "max_spread_direction": u_S,
        }
    else:
        return {
            "phi_E"               : phi_E,
            "max_spread_direction": np.asarray((0,1,0)), # default: North
        }
# surface-fire-combine-wind-and-slope-vectors ends here
# [[file:../../org/pyretechnics.org::surface-fire-eccentricity][surface-fire-eccentricity]]
from math import exp, sqrt
from pyretechnics.conversion import m_min_to_mph


def surface_length_to_width_ratio(effective_wind_speed, model="rothermel"):
    """
    Calculate the length_to_width_ratio of the surface fire front given:
    - effective_wind_speed :: m/min (aligned with the slope-tangential plane)
    - model                :: "rothermel" or "behave" (Optional)
    """
    effective_wind_speed_mph = m_min_to_mph(effective_wind_speed)
    match model:
        case "rothermel":
            return 1.0 + 0.25 * effective_wind_speed_mph

        case "behave":
            return min(8.0,
                       0.936 * exp(0.1147 * effective_wind_speed_mph)
                       +
                       0.461 * exp(-0.0692 * effective_wind_speed_mph)
                       -
                       0.397)

        case _:
            raise ValueError("Invalid input: model must be 'rothermel' or 'behave'.")


def surface_fire_eccentricity(length_to_width_ratio):
    """
    Calculate the eccentricity (E) of the surface fire front using eq. 8 from
    Albini and Chase 1980 given:
    - L/W :: (1: circular spread, > 1: elliptical spread)
    """
    return sqrt(length_to_width_ratio ** 2.0 - 1.0) / length_to_width_ratio
# surface-fire-eccentricity ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-max][surface-fire-behavior-max]]
from pyretechnics.conversion import opposite_direction
from pyretechnics.vector_utils import vector_magnitude


def maybe_limit_wind_speed(use_wind_limit, max_wind_speed, get_phi_W, get_wind_speed, phi_E_magnitude):
    """
    Given these inputs:
    - use_wind_limit  :: boolean
    - max_wind_speed  :: m/min
    - get_phi_W       :: lambda: midflame_wind_speed (m/min) => phi_W (unitless)
    - get_wind_speed  :: lambda: phi_W (unitless) => midflame_wind_speed (m/min)
    - phi_E_magnitude :: unitless

    return a tuple with these fields:
    - limited_wind_speed :: m/min
    - limited_phi_E      :: unitless
    """
    effective_wind_speed = get_wind_speed(phi_E_magnitude)
    if (use_wind_limit and effective_wind_speed > max_wind_speed):
        return (
            max_wind_speed,
            get_phi_W(max_wind_speed),
        )
    else:
        return (
            effective_wind_speed,
            phi_E_magnitude,
        )


# NOTE: No longer takes ellipse_adjustment_factor parameter
def calc_surface_fire_behavior_max(surface_fire_min, midflame_wind_speed, upwind_direction,
                                   slope, aspect, use_wind_limit=True):
    """
    Given these inputs:
    - surface_fire_min       :: dictionary of no-wind-no-slope surface fire behavior values
      - base_spread_rate         :: m/min
      - base_fireline_intensity  :: kW/m
      - max_effective_wind_speed :: m/min
      - get_phi_S                :: lambda: slope (rise/run) => phi_S (unitless)
      - get_phi_W                :: lambda: midflame_wind_speed (m/min) => phi_W (unitless)
      - get_wind_speed           :: lambda: phi_W (unitless) => midflame_wind_speed (m/min)
    - midflame_wind_speed    :: m/min
    - upwind_direction       :: degrees clockwise from North
    - slope                  :: rise/run
    - aspect                 :: degrees clockwise from North
    - use_wind_limit         :: boolean (Optional)

    return a dictionary containing these keys:
    - max_spread_rate        :: m/min
    - max_spread_direction   :: (x, y, z) unit vector
    - max_fireline_intensity :: kW/m
    - max_flame_length       :: m
    - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
    - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
    """
    # Unpack no-wind-no-slope surface fire behavior values
    spread_rate        = surface_fire_min["base_spread_rate"]
    fireline_intensity = surface_fire_min["base_fireline_intensity"]
    max_wind_speed     = surface_fire_min["max_effective_wind_speed"]
    get_phi_W          = surface_fire_min["get_phi_W"]
    get_phi_S          = surface_fire_min["get_phi_S"]
    get_wind_speed     = surface_fire_min["get_wind_speed"]
    # Reverse the provided wind and slope directions
    downwind_direction = opposite_direction(upwind_direction)
    upslope_direction  = opposite_direction(aspect)
    # Project wind and slope vectors onto the slope-tangential plane
    vectors = project_wind_and_slope_vectors_3d(midflame_wind_speed, downwind_direction, slope, upslope_direction)
    wind_vector_3d  = vectors["wind_vector_3d"]  # m/min
    slope_vector_3d = vectors["slope_vector_3d"] # rise/run
    # Calculate phi_W and phi_S
    phi_W = get_phi_W(vector_magnitude(wind_vector_3d)) # |wind_vector_3d| = slope-aligned midflame wind speed
    phi_S = get_phi_S(slope)
    # Calculate phi_E and the max_spread_direction
    result               = get_phi_E(wind_vector_3d, slope_vector_3d, phi_W, phi_S)
    phi_E                = result["phi_E"]
    max_spread_direction = result["max_spread_direction"]
    # Limit effective wind speed to max wind speed if use_wind_limit == True
    (limited_wind_speed, limited_phi_E) = maybe_limit_wind_speed(use_wind_limit, max_wind_speed,
                                                                 get_phi_W, get_wind_speed, phi_E)
    # Calculate and return max surface fire behavior values
    max_spread_rate        = spread_rate * (1.0 + limited_phi_E)
    max_fireline_intensity = fireline_intensity * (1.0 + limited_phi_E)
    length_to_width_ratio  = surface_length_to_width_ratio(limited_wind_speed)
    return {
        "max_spread_rate"       : max_spread_rate,
        "max_spread_direction"  : max_spread_direction, # unit vector
        "max_fireline_intensity": max_fireline_intensity,
        "max_flame_length"      : calc_flame_length(max_fireline_intensity),
        "length_to_width_ratio" : length_to_width_ratio,
        "eccentricity"          : surface_fire_eccentricity(length_to_width_ratio),
    }
# surface-fire-behavior-max ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-in-direction][surface-fire-behavior-in-direction]]
import numpy as np


def calc_surface_fire_behavior_in_direction(surface_fire_max, spread_direction):
    """
    Given these inputs:
    - surface_fire_max     :: dictionary of max surface fire behavior values
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_fireline_intensity :: kW/m
      - max_flame_length       :: m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
    - spread_direction     :: 3D unit vector on the slope-tangential plane

    return a dictionary containing these keys:
    - fire_type          :: "surface"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    # Unpack max surface fire behavior values
    max_spread_rate        = surface_fire_max["max_spread_rate"]
    max_spread_direction   = surface_fire_max["max_spread_direction"]
    max_fireline_intensity = surface_fire_max["max_fireline_intensity"]
    eccentricity           = surface_fire_max["eccentricity"]
    # Calculate cos(w), where w is the offset angle between these unit vectors on the slope-tangential plane
    cos_w = np.dot(max_spread_direction, np.asarray(spread_direction))
    # Calculate adjustment due to the offset angle from the max spread direction
    adjustment = (1.0 - eccentricity) / (1.0 - eccentricity * cos_w)
    # Update surface fire behavior values by the adjustment value
    fireline_intensity = max_fireline_intensity * adjustment
    return {
        "fire_type"         : "surface",
        "spread_rate"       : max_spread_rate * adjustment,
        "spread_direction"  : spread_direction,
        "fireline_intensity": fireline_intensity,
        "flame_length"      : calc_flame_length(fireline_intensity),
    }
# surface-fire-behavior-in-direction ends here
