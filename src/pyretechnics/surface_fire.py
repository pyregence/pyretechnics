# [[file:../../org/Pyretechnics.org::rothermel-surface-fire-spread-no-wind-no-slope][rothermel-surface-fire-spread-no-wind-no-slope]]
from math import exp
from fuel_models import map_category, map_size_class, category_sum, size_class_sum

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


def calc_packing_ratio(w_o, rho_p, delta):
    if (delta > 0.0):
        beta_i = size_class_sum(lambda i: w_o[i] / rho_p[i])
        return category_sum(lambda i: beta_i[i]) / delta
    else:
        return 0.0


def calc_surface_area_to_volume_ratio(f_i, f_ij, sigma):
    sigma_prime_i = size_class_sum(lambda i: f_ij[i] * sigma[i])
    return category_sum(lambda i: f_i[i] * sigma_prime_i[i])


def calc_optimum_packing_ratio(sigma_prime):
    return (3.348 / sigma_prime ** 0.8189) if (sigma_prime > 0.0) else 1.0


def calc_optimum_reaction_velocity(beta, sigma_prime, beta_op):
    # Albini 1976 replaces 1 / (4.774 * (sigma_prime ** 0.1) - 7.27)
    A               = (133.0 / sigma_prime ** 0.7913) if (sigma_prime > 0.0) else 0.0
    B               = sigma_prime ** 1.5
    C               = beta / beta_op
    # Maximum reaction velocity (1/min)
    Gamma_prime_max = B / (495.0 + 0.0594 * B)
    # Optimum reaction velocity (1/min)
    return Gamma_prime_max * (C ** A) * exp(A * (1.0 - C))


def calc_heat_per_unit_area(eta_S_i, eta_M_i, h_i, W_n_i):
    return category_sum(lambda i: W_n_i[i] * h_i[i] * eta_M_i[i] * eta_S_i[i])


def calc_reaction_intensity(Gamma_prime, Btus):
    return Gamma_prime * Btus


def calc_propagating_flux_ratio(beta, sigma_prime):
    return exp((0.792 + 0.681 * (sigma_prime ** 0.5)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma_prime)


def calc_heat_of_preignition(M_f):
    return map_size_class(lambda i: 250.0 + 1116.0 * M_f[i])


def calc_heat_distribution(sigma, Q_ig, f_ij):
    return size_class_sum(lambda i:
                          (lambda sigma:
                           f_ij[i] * exp(-138.0 / sigma) * Q_ig[i] if (sigma > 0.0) else 0.0
                           )(sigma[i]))


def calc_ovendry_bulk_density(w_o, delta):
    if (delta > 0.0):
        rho_b_i = size_class_sum(lambda i: w_o[i])
        return category_sum(lambda i: rho_b_i[i]) / delta
    else:
        return 0.0


def calc_heat_total(f_i, epsilon_i):
    return category_sum(lambda i: f_i[i] * epsilon_i[i])


def calc_surface_fire_spread_rate(I_R, xi, rho_b, epsilon):
    rho_b_epsilon_Q_ig = rho_b * epsilon
    return ((I_R * xi) / rho_b_epsilon_Q_ig) if (rho_b_epsilon_Q_ig > 0.0) else 0.0


def calc_residence_time(sigma_prime):
    return 384.0 / sigma_prime


def get_phi_S_fn(beta):
    if (beta > 0.0):
        G = 5.275 * beta ** -0.3
        return lambda slope: (slope ** 2.0) * G if (slope > 0.0) else 0.0
    else:
        return lambda _: 0.0


def get_phi_W_fn(beta, B, C, F):
    if (beta > 0.0):
        C_over_F = C / F
        return lambda midflame_wind_speed: (midflame_wind_speed ** B) * C_over_F if (midflame_wind_speed > 0.0) else 0.0
    else:
        return lambda _: 0.0


def get_wind_speed_fn(B, C, F):
    F_over_C  = F / C
    B_inverse = 1.0 / B
    return lambda phi_W: (phi_W * F_over_C) ** B_inverse


def rothermel_surface_fire_spread_no_wind_no_slope(fuel_model):
    """
    Returns the rate of surface fire spread in ft/min and the reaction
    intensity (i.e., amount of heat output) of a fire in Btu/ft^2*min
    given a map containing these keys:
    - delta [fuel depth (ft)]
    - w_o [ovendry fuel loading (lb/ft^2)]
    - sigma [fuel particle surface-area-to-volume ratio (ft^2/ft^3)]
    - h [fuel particle low heat content (Btu/lb)]
    - rho_p [ovendry particle density (lb/ft^3)]
    - S_T [fuel particle total mineral content (lb minerals/lb ovendry wood)]
    - S_e [fuel particle effective mineral content (lb silica-free minerals/lb ovendry wood)]
    - M_x [moisture content of extinction (lb moisture/lb ovendry wood)]
    - M_f [fuel particle moisture content (lb moisture/lb ovendry wood)]
    - f_ij [percent of load per size class (%)]
    - f_i [percent of load per category (%)]
    - g_ij [percent of load per size class from Albini_1976_FIREMOD, page 20]
    """
    delta          = fuel_model["delta"]
    w_o            = fuel_model["w_o"]
    sigma          = fuel_model["sigma"]
    h              = fuel_model["h"]
    rho_p          = fuel_model["rho_p"]
    S_T            = fuel_model["S_T"]
    S_e            = fuel_model["S_e"]
    M_x            = fuel_model["M_x"]
    M_f            = fuel_model["M_f"]
    f_ij           = fuel_model["f_ij"]
    f_i            = fuel_model["f_i"]
    g_ij           = fuel_model["g_ij"]
    eta_S_i        = calc_mineral_damping_coefficients(f_ij, S_e)
    eta_M_i        = calc_moisture_damping_coefficients(f_ij, M_f, M_x)
    h_i            = calc_low_heat_content(f_ij, h)
    W_n_i          = calc_net_fuel_loading(g_ij, w_o, S_T)                      # (lb/ft^2)
    beta           = calc_packing_ratio(w_o, rho_p, delta)
    sigma_prime    = calc_surface_area_to_volume_ratio(f_i, f_ij, sigma)
    beta_op        = calc_optimum_packing_ratio(sigma_prime)
    Gamma_prime    = calc_optimum_reaction_velocity(beta, sigma_prime, beta_op) # (1/min)
    Btus           = calc_heat_per_unit_area(eta_S_i, eta_M_i, h_i, W_n_i)      # (Btu/ft^2)
    I_R            = calc_reaction_intensity(Gamma_prime, Btus)                 # (Btu/ft^2*min)
    xi             = calc_propagating_flux_ratio(beta, sigma_prime)
    Q_ig           = calc_heat_of_preignition(M_f)                              # (Btu/lb)
    epsilon_i      = calc_heat_distribution(sigma, Q_ig, f_ij)                  # (Btu/lb)
    rho_b          = calc_ovendry_bulk_density(w_o, delta)                      # (lb/ft^3)
    epsilon        = calc_heat_total(f_i, epsilon_i)                            # (Btu/lb)
    R              = calc_surface_fire_spread_rate(I_R, xi, rho_b, epsilon)     # (ft/min)
    t_res          = calc_residence_time(sigma_prime)
    B              = 0.02526 * (sigma_prime ** 0.54)
    C              = 7.47 * exp(-0.133 * (sigma_prime ** 0.55))
    E              = 0.715 * exp(-3.59 * (sigma_prime / 10000.0))
    F              = (beta / beta_op) ** E
    get_phi_S      = get_phi_S_fn(beta)
    get_phi_W      = get_phi_W_fn(beta, B, C, F)
    get_wind_speed = get_wind_speed_fn(B, C, F)
    return {
        "unadj_spread_rate" : R,
        "reaction_intensity": I_R,
        "residence_time"    : t_res,
        "fuel_bed_depth"    : delta,
        "heat_of_combustion": h[0],
        "get_phi_S"         : get_phi_S,
        "get_phi_W"         : get_phi_W,
        "get_wind_speed"    : get_wind_speed,
    }

# Test with:
# - moisturize(fuel_models_precomputed[1], [0.05, 0.10, 0.15, 0.05, 0.30, 0.50])   # non-dynamic fuel model
# - moisturize(fuel_models_precomputed[101], [0.05, 0.10, 0.15, 0.05, 0.30, 0.50]) # dynamic fuel model
# rothermel-surface-fire-spread-no-wind-no-slope ends here
# [[file:../../org/Pyretechnics.org::wind-adjustment-factor][wind-adjustment-factor]]
from math import log, sqrt

def wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover):
    """
    fuel_bed_depth :: ft
    canopy_height  :: ft
    canopy_cover   :: 0-100
    """
    if (canopy_cover > 0.0) and (canopy_height > 0.0):
        # sheltered: equation 2 based on CC and CH, CR=1 (Andrews 2012)
        A = sqrt((canopy_cover / 300.0) * canopy_height)
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
