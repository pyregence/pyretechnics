# [[file:../../org/pyretechnics.org::surface-fire-imports][surface-fire-imports]]
import cython
import cython as cy
if cython.compiled:
    from cython.cimports.libc.math import sqrt, pow, log, exp
    from cython.cimports.pyretechnics.cy_types import \
        vec_xy, vec_xyz, fcatarr, fclaarr, FuelModel, ProjectedVectors, FireBehaviorMin, \
        FireBehaviorMax, SpreadBehavior
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.vector_utils as vu
else:
    from math import sqrt, pow, log, exp
    from pyretechnics.py_types import \
        vec_xy, vec_xyz, fcatarr, fclaarr, FuelModel, ProjectedVectors, FireBehaviorMin, \
        FireBehaviorMax, SpreadBehavior
    import pyretechnics.conversion as conv
    import pyretechnics.vector_utils as vu
# surface-fire-imports ends here
# [[file:../../org/pyretechnics.org::surface-fire-utility-functions][surface-fire-utility-functions]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def __dotp_in_category(x_ij: fclaarr, y_ij: fclaarr) -> fcatarr:
    """
    Util: dot-product between fuel size class vectors, within each fuel category (dead/live).
    """
    return (
        (
            (x_ij[0] * y_ij[0]) +
            (x_ij[1] * y_ij[1]) +
            (x_ij[2] * y_ij[2]) +
            (x_ij[3] * y_ij[3])
        ),
        (
            (x_ij[4] * y_ij[4]) +
            (x_ij[5] * y_ij[5])
        )
    )


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def __dotp_categories(x_i: fcatarr, y_i: fcatarr) -> cy.float:
    return (
        (x_i[0] * y_i[0]) +
        (x_i[1] * y_i[1])
    )
# surface-fire-utility-functions ends here
# [[file:../../org/pyretechnics.org::surface-fire-common-intermediate-calculations][surface-fire-common-intermediate-calculations]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_surface_area_to_volume_ratio(f_i: fcatarr, f_ij: fclaarr, sigma: fclaarr) -> cy.float:
    sigma_prime: fcatarr = __dotp_in_category(f_ij, sigma)
    return __dotp_categories(f_i, sigma_prime)


@cy.cfunc
@cy.exceptval(check=False)
def calc_packing_ratio(w_o: fclaarr, rho_p: fclaarr, delta: cy.float) -> cy.float:
    if (delta > 0.0):
        rho_p_inv: fclaarr = ( # TODO: OPTIM pre-compute
            1.0/rho_p[0],
            1.0/rho_p[1],
            1.0/rho_p[2],
            1.0/rho_p[3],
            1.0/rho_p[4],
            1.0/rho_p[5]
        )
        beta: fcatarr = __dotp_in_category(w_o, rho_p_inv)
        return (beta[0] + beta[1]) / delta
    else:
        return 0.0


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_optimum_packing_ratio(sigma_prime: cy.float) -> cy.float:
    return (3.348 * pow(sigma_prime, -0.8189)) if (sigma_prime > 0.0) else 1.0
# surface-fire-common-intermediate-calculations ends here
# [[file:../../org/pyretechnics.org::surface-fire-reaction-intensity][surface-fire-reaction-intensity]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_mineral_damping_coefficients(f_ij: fclaarr, S_e: fclaarr) -> fcatarr:
    S_e_i: fcatarr = __dotp_in_category(f_ij, S_e)
    (S_e_0, S_e_1) = S_e_i
    return (
        0.174 * pow(S_e_0, -0.19) if (S_e_0 > 0.0) else 1.0,
        0.174 * pow(S_e_1, -0.19) if (S_e_1 > 0.0) else 1.0
    )


@cy.cfunc
@cy.exceptval(check=False)
def __cat_moisture_damping_coefficient(M_f: cy.float, M_x: cy.float) -> cy.float:
    if (M_x > 0.0):
        r_M : cy.float = min(1.0, M_f / M_x)
        r_M2: cy.float = r_M * r_M
        r_M3: cy.float = r_M2 * r_M
        return 1.0 - (2.59 * r_M) + (5.11 * r_M2) - (3.52 * r_M3)
    else:
        return 0.0


@cy.cfunc
@cy.exceptval(check=False)
def calc_moisture_damping_coefficients(f_ij: fclaarr, M_f: fclaarr, M_x: fclaarr) -> fcatarr:
    M_f_i: fcatarr = __dotp_in_category(f_ij, M_f)
    M_x_i: fcatarr = __dotp_in_category(f_ij, M_x)
    return (
        __cat_moisture_damping_coefficient(M_f_i[0], M_x_i[0]),
        __cat_moisture_damping_coefficient(M_f_i[1], M_x_i[1])
    )


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_low_heat_content(f_ij: fclaarr, h: fclaarr) -> fcatarr:
    return __dotp_in_category(f_ij, h)


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_net_fuel_loading(g_ij: fclaarr, w_o: fclaarr, S_T: fclaarr) -> fcatarr:
    return (
        (
            g_ij[0] * w_o[0] * (1.0 - S_T[0]) +
            g_ij[1] * w_o[1] * (1.0 - S_T[1]) +
            g_ij[2] * w_o[2] * (1.0 - S_T[2]) +
            g_ij[3] * w_o[3] * (1.0 - S_T[3])
        ),
        (
            g_ij[4] * w_o[4] * (1.0 - S_T[4]) +
            g_ij[5] * w_o[5] * (1.0 - S_T[5])
        )
    )


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_heat_per_unit_area(eta_S_i: fcatarr, eta_M_i: fcatarr, h_i: fcatarr, W_n_i: fcatarr) -> cy.float:
    return (
        (W_n_i[0] * h_i[0] * eta_M_i[0] * eta_S_i[0]) +
        (W_n_i[1] * h_i[1] * eta_M_i[1] * eta_S_i[1])
    )


@cy.cfunc
@cy.exceptval(check=False)
def calc_optimum_reaction_velocity(sigma_prime: cy.float, beta: cy.float, beta_op: cy.float) -> cy.float:
    # Albini 1976 replaces 1 / (4.774 * (sigma_prime ** 0.1) - 7.27)
    A: cy.float = (133.0 * pow(sigma_prime, -0.7913)) if (sigma_prime > 0.0) else 0.0
    B: cy.float = pow(sigma_prime, 1.5)
    C: cy.float = beta / beta_op
    # Maximum reaction velocity (1/min)
    Gamma_prime_max: cy.float = B / (495.0 + 0.0594 * B)
    # Optimum reaction velocity (1/min)
    return Gamma_prime_max * pow(C, A) * exp(A * (1.0 - C))


@cy.cfunc
@cy.exceptval(check=False)
def calc_reaction_intensity(moisturized_fuel_model: FuelModel,
                            sigma_prime           : cy.float,
                            beta                  : cy.float,
                            beta_op               : cy.float) -> cy.float:
    # Unpack the FuelModel fields
    w_o : fclaarr = moisturized_fuel_model.w_o
    h   : fclaarr = moisturized_fuel_model.h
    S_T : fclaarr = moisturized_fuel_model.S_T
    S_e : fclaarr = moisturized_fuel_model.S_e
    M_x : fclaarr = moisturized_fuel_model.M_x
    M_f : fclaarr = moisturized_fuel_model.M_f
    f_ij: fclaarr = moisturized_fuel_model.f_ij
    g_ij: fclaarr = moisturized_fuel_model.g_ij
    # Derive intermediate quantities
    eta_S_i    : fcatarr  = calc_mineral_damping_coefficients(f_ij, S_e)
    eta_M_i    : fcatarr  = calc_moisture_damping_coefficients(f_ij, M_f, M_x)
    h_i        : fcatarr  = calc_low_heat_content(f_ij, h)                             # (Btu/lb)
    W_n_i      : fcatarr  = calc_net_fuel_loading(g_ij, w_o, S_T)                      # (lb/ft^2)
    Btus       : cy.float = calc_heat_per_unit_area(eta_S_i, eta_M_i, h_i, W_n_i)      # (Btu/ft^2)
    Gamma_prime: cy.float = calc_optimum_reaction_velocity(sigma_prime, beta, beta_op) # (1/min)
    # Calculate reaction intensity
    return Btus * Gamma_prime # (Btu/ft^2/min)
# surface-fire-reaction-intensity ends here
# [[file:../../org/pyretechnics.org::surface-fire-propagating-flux-ratio][surface-fire-propagating-flux-ratio]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_propagating_flux_ratio(sigma_prime: cy.float, beta: cy.float) -> cy.float:
    return exp((0.792 + 0.681 * sqrt(sigma_prime)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma_prime)
# surface-fire-propagating-flux-ratio ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-source-no-wind-no-slope][surface-fire-heat-source-no-wind-no-slope]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_heat_source(I_R: cy.float, xi: cy.float) -> cy.float:
    return I_R * xi
# surface-fire-heat-source-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::surface-fire-oven-dry-fuel-bed-bulk-density][surface-fire-oven-dry-fuel-bed-bulk-density]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_ovendry_bulk_density(w_o: fclaarr, delta: cy.float) -> cy.float:
    if (delta > 0.0):
        w_o_sum: cy.float = w_o[0] + w_o[1] + w_o[2] + w_o[3] + w_o[4] + w_o[5]
        return w_o_sum / delta
    else:
        return 0.0
# surface-fire-oven-dry-fuel-bed-bulk-density ends here
# [[file:../../org/pyretechnics.org::surface-fire-effective-heating-number-distribution][surface-fire-effective-heating-number-distribution]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def __sizeclass_heating_number(sigma_i: cy.float) -> cy.float:
    return exp(-138.0 / sigma_i) if (sigma_i > 0.0) else 0.0


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_effective_heating_number_distribution(sigma: fclaarr) -> fclaarr: # OPTIM pre-compute, exp is expensive
    return (
        __sizeclass_heating_number(sigma[0]),
        __sizeclass_heating_number(sigma[1]),
        __sizeclass_heating_number(sigma[2]),
        __sizeclass_heating_number(sigma[3]),
        __sizeclass_heating_number(sigma[4]),
        __sizeclass_heating_number(sigma[5])
    )
# surface-fire-effective-heating-number-distribution ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-of-preignition-distribution][surface-fire-heat-of-preignition-distribution]]
@cy.cfunc
@cy.exceptval(check=False)
def __sizeclass_heat_of_preignition_distribution(M_f_i: cy.float) -> cy.float:
    pos1: cy.float = cy.cast(cy.float, (M_f_i > 0.0))
    return (250.0 + 1116.0 * M_f_i) * pos1 # Returns 0 unless M_f_i > 0


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_heat_of_preignition_distribution(M_f: fclaarr) -> fclaarr:
    return (
        __sizeclass_heat_of_preignition_distribution(M_f[0]),
        __sizeclass_heat_of_preignition_distribution(M_f[1]),
        __sizeclass_heat_of_preignition_distribution(M_f[2]),
        __sizeclass_heat_of_preignition_distribution(M_f[3]),
        __sizeclass_heat_of_preignition_distribution(M_f[4]),
        __sizeclass_heat_of_preignition_distribution(M_f[5])
    )
# surface-fire-heat-of-preignition-distribution ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-sink][surface-fire-heat-sink]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_heat_sink(f_i: fcatarr, f_ij: fclaarr, rho_b: cy.float, epsilon_ij: fclaarr, Q_ig_ij: fclaarr) -> cy.float:
    effective_heat_of_preignition_i: fcatarr = (
        (
            (f_ij[0] * epsilon_ij[0] * Q_ig_ij[0]) +
            (f_ij[1] * epsilon_ij[1] * Q_ig_ij[1]) +
            (f_ij[2] * epsilon_ij[2] * Q_ig_ij[2]) +
            (f_ij[3] * epsilon_ij[3] * Q_ig_ij[3])
        ),
        (
            (f_ij[4] * epsilon_ij[4] * Q_ig_ij[4]) +
            (f_ij[5] * epsilon_ij[5] * Q_ig_ij[5])
        )
    )
    effective_heat_of_preignition: cy.float = __dotp_categories(f_i, effective_heat_of_preignition_i)
    return rho_b * effective_heat_of_preignition
# surface-fire-heat-sink ends here
# [[file:../../org/pyretechnics.org::surface-fire-spread-rate-no-wind-no-slope][surface-fire-spread-rate-no-wind-no-slope]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_spread_rate(heat_source: cy.float, heat_sink: cy.float) -> cy.float:
    return heat_source / heat_sink if (heat_sink > 0.0) else 0.0
# surface-fire-spread-rate-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::surface-fire-intensity-functions][surface-fire-intensity-functions]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_residence_time(sigma_prime: cy.float) -> cy.float:
    """
    Returns the residence time (total burning time) of fuel (min) given:
    - sigma_prime :: ft^2/ft^3 (surface area to volume ratio)
    """
    return 384.0 / sigma_prime if (sigma_prime > 0.0) else 0.0


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_flame_depth(spread_rate: cy.float, residence_time: cy.float) -> cy.float:
    """
    Returns the depth, or front-to-back distance, of the actively flaming zone
    of a free-spreading fire (ft) given:
    - spread_rate    :: ft/min (orthogonal to the fireline)
    - residence_time :: min
    """
    return spread_rate * residence_time


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_fireline_intensity(reaction_intensity: cy.float, flame_depth: cy.float) -> cy.float:
    """
    Returns the rate of heat release per unit of fire edge (Btu/ft/s) given:
    - reaction_intensity :: Btu/ft^2/min
    - flame_depth        :: ft
    """
    return (reaction_intensity * flame_depth) / 60.0


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_flame_length(fireline_intensity: cy.float) -> cy.float:
    """
    Returns the average flame length (m) given:
    - fireline_intensity :: kW/m
    """
    return 0.07747042253266703 * pow(fireline_intensity, 0.46)


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def calc_areal_heat_output(spread_rate: cy.float, fireline_intensity: cy.float) -> cy.float:
    """
    Returns the heat per unit area (kJ/m^2) given:
    - spread_rate        :: m/min
    - fireline_intensity :: kW/m
    """
    return 60.0 * fireline_intensity / spread_rate if spread_rate > 0.0 else 0.0
# surface-fire-intensity-functions ends here
# [[file:../../org/pyretechnics.org::surface-fire-max-effective-wind-speed][surface-fire-max-effective-wind-speed]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_max_effective_wind_speed(reaction_intensity: cy.float) -> cy.float:
    return 0.9 * reaction_intensity
# surface-fire-max-effective-wind-speed ends here
# [[file:../../org/pyretechnics.org::surface-fire-slope-factor-function][surface-fire-slope-factor-function]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def get_phi_S(surface_fire_min: FireBehaviorMin, slope: cy.float) -> cy.float:
    return surface_fire_min._phiS_G * (slope * slope)
# surface-fire-slope-factor-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-wind-factor-function][surface-fire-wind-factor-function]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def get_phi_W(surface_fire_min: FireBehaviorMin, midflame_wind_speed: cy.float) -> cy.float:
    return surface_fire_min._phiW_scalr * pow(midflame_wind_speed, surface_fire_min._phiW_expnt)
# surface-fire-wind-factor-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-wind-speed-function][surface-fire-wind-speed-function]]
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def get_wind_speed(surface_fire_min: FireBehaviorMin, phi_W: cy.float) -> cy.float:
    return surface_fire_min._ws_scalr * pow(phi_W, surface_fire_min._ws_expnt)
# surface-fire-wind-speed-function ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-no-wind-no-slope][surface-fire-behavior-no-wind-no-slope]]
@cy.cfunc
@cy.exceptval(check=False)
def make_surface_fire_min(base_spread_rate        : cy.float,
                          base_fireline_intensity : cy.float,
                          max_effective_wind_speed: cy.float,
                          B                       : cy.float,
                          C                       : cy.float,
                          F                       : cy.float,
                          beta                    : cy.float) -> FireBehaviorMin:
    _phiS_G: cy.float = 0.0
    if (beta > 0.0):
        _phiS_G = 5.275 * pow(beta, -0.3)

    _phiW_scalr: cy.float = 0.0
    _phiW_expnt: cy.float = 0.0
    if (F > 0.0):
        _phiW_scalr = (C / F) * pow(conv.m_to_ft(1.0), B)
        _phiW_expnt = B

    _ws_scalr: cy.float = 0.0
    _ws_expnt: cy.float = 0.0
    if (B > 0.0):
        B_inverse: cy.float = 1.0 / B
        _ws_scalr = conv.ft_to_m(pow((F / C), B_inverse))
        _ws_expnt = B_inverse

    return FireBehaviorMin(
        base_spread_rate,
        base_fireline_intensity,
        max_effective_wind_speed,
        _phiS_G,
        _phiW_scalr,
        _phiW_expnt,
        _ws_scalr,
        _ws_expnt)


@cy.ccall
@cy.exceptval(check=False)
def calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model: FuelModel,
                                                spread_rate_adjustment: cy.float = 1.0) -> FireBehaviorMin:
    """
    Given these inputs:
    - moisturized_fuel_model :: FuelModel struct of fuel model and fuel moisture properties
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

    return a FireBehaviorMin struct containing these keys:
    - base_spread_rate         :: m/min
    - base_fireline_intensity  :: kW/m
    - max_effective_wind_speed :: m/min
    - _phiS_G                  :: intermediate value for computing phi_S (unitless)
    - _phiW_scalr              :: intermediate value for computing phi_W (unitless)
    - _phiW_expnt              :: intermediate value for computing phi_W (unitless)
    - _ws_scalr                :: intermediate value for computing effective_wind_speed (m/min)
    - _ws_expnt                :: intermediate value for computing effective_wind_speed (unitless)
    """
    # Unpack fuel model values
    delta: cy.float = moisturized_fuel_model.delta
    w_o  : fclaarr  = moisturized_fuel_model.w_o
    rho_p: fclaarr  = moisturized_fuel_model.rho_p
    sigma: fclaarr  = moisturized_fuel_model.sigma
    M_f  : fclaarr  = moisturized_fuel_model.M_f
    f_ij : fclaarr  = moisturized_fuel_model.f_ij
    f_i  : fcatarr  = moisturized_fuel_model.f_i
    # Calculate base spread rate (no wind, no slope)
    sigma_prime: cy.float = calc_surface_area_to_volume_ratio(f_i, f_ij, sigma)
    beta       : cy.float = calc_packing_ratio(w_o, rho_p, delta)
    beta_op    : cy.float = calc_optimum_packing_ratio(sigma_prime)
    I_R        : cy.float = calc_reaction_intensity(moisturized_fuel_model, sigma_prime, beta, beta_op) # Btu/ft^2/min
    xi         : cy.float = calc_propagating_flux_ratio(sigma_prime, beta)
    heat_source: cy.float = calc_heat_source(I_R, xi)                             # Btu/ft^2/min
    rho_b      : cy.float = calc_ovendry_bulk_density(w_o, delta)                 # lb/ft^3
    epsilon_ij : fclaarr  = calc_effective_heating_number_distribution(sigma)
    Q_ig_ij    : fclaarr  = calc_heat_of_preignition_distribution(M_f)            # Btu/lb
    heat_sink  : cy.float = calc_heat_sink(f_i, f_ij, rho_b, epsilon_ij, Q_ig_ij) # Btu/ft^3
    R0         : cy.float = calc_spread_rate(heat_source, heat_sink)              # ft/min
    # Calculate base fireline intensity (no wind, no slope)
    t_res: cy.float = calc_residence_time(sigma_prime)  # min
    D_A  : cy.float = calc_flame_depth(R0, t_res)       # ft
    I_s  : cy.float = calc_fireline_intensity(I_R, D_A) # Btu/ft/s
    # Pre-compute values related to wind and slope
    U_eff_max: cy.float = calc_max_effective_wind_speed(I_R) # ft/min
    B        : cy.float = 0.02526 * pow(sigma_prime, 0.54)
    C        : cy.float = 7.47 * exp(-0.133 * pow(sigma_prime, 0.55))
    E        : cy.float = 0.715 * exp(-3.59 * (sigma_prime * 1e-4))
    F        : cy.float = pow((beta / beta_op), E)
    # Return no-wind-no-slope surface fire behavior values
    base_spread_rate        : cy.float = conv.ft_to_m(R0 * spread_rate_adjustment)
    base_fireline_intensity : cy.float = conv.Btu_ft_s_to_kW_m(I_s * spread_rate_adjustment)
    max_effective_wind_speed: cy.float = conv.ft_to_m(U_eff_max)
    return make_surface_fire_min(
        base_spread_rate,
        base_fireline_intensity,
        max_effective_wind_speed,
        B,
        C,
        F,
        beta
    )
# surface-fire-behavior-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::midflame-wind-speed][midflame-wind-speed]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_wind_adjustment_factor(fuel_bed_depth: cy.float, canopy_height: cy.float, canopy_cover: cy.float) -> cy.float:
    """
    Return the wind adjustment factor (unitless) given these inputs:
    - fuel_bed_depth :: ft
    - canopy_height  :: ft
    - canopy_cover   :: 0-1
    """
    if (canopy_cover > 0.0) and (canopy_height > 0.0):
        # sheltered: equation 2 based on CC and CH, CR=1 (Andrews 2012)
        A: cy.float = sqrt((canopy_cover / 3.0) * canopy_height)
        B: cy.float = log((20.0 + 0.36 * canopy_height) / (0.13 * canopy_height))
        return 0.555 / (A * B)
    elif (fuel_bed_depth > 0.0):
        # unsheltered: equation 6 H_F = H (Andrews 2012)
        A: cy.float = log((20.0 + 0.36 * fuel_bed_depth) / (0.13 * fuel_bed_depth))
        return 1.83 / A # 1.83 truncated from 1.8328795184533409
    else:
        # non-burnable fuel model
        return 0.0


@cy.ccall
@cy.exceptval(check=False)
def calc_midflame_wind_speed(wind_speed_20ft: cy.float,
                             fuel_bed_depth : cy.float,
                             canopy_height  : cy.float,
                             canopy_cover   : cy.float) -> cy.float:
    """
    Return the midflame wind speed (S) given these inputs:
    - wind_speed_20ft :: S
    - fuel_bed_depth  :: ft
    - canopy_height   :: ft
    - canopy_cover    :: 0-1
    """
    wind_adj_factor: cy.float = calc_wind_adjustment_factor(fuel_bed_depth, canopy_height, canopy_cover)
    return wind_speed_20ft * wind_adj_factor
# midflame-wind-speed ends here
# [[file:../../org/pyretechnics.org::surface-fire-combine-wind-and-slope-vectors][surface-fire-combine-wind-and-slope-vectors]]
@cy.cfunc
@cy.exceptval(check=False)
def project_wind_and_slope_vectors_3d(wind_speed        : cy.float,
                                      downwind_direction: cy.float,
                                      slope             : cy.float,
                                      upslope_direction : cy.float) -> ProjectedVectors:
    """
    Given these inputs:
    - wind_speed         :: S
    - downwind_direction :: degrees clockwise from North
    - slope              :: rise/run
    - upslope_direction  :: degrees clockwise from North

    return a ProjectedVectors struct containing these keys:
    - wind_vector_3d  :: (x: S, y: S, z: S)
    - slope_vector_3d :: (x, y, z)
    """
    # Convert wind and slope vectors from azimuthal to cartesian coordinates
    wind_vector_2d : vec_xy = conv.azimuthal_to_cartesian(wind_speed, downwind_direction)
    slope_vector_2d: vec_xy = conv.azimuthal_to_cartesian(slope, upslope_direction)
    # Project wind and slope vectors onto the slope-tangential plane
    wind_vector_3d : vec_xyz = vu.to_slope_plane(wind_vector_2d, slope_vector_2d)
    slope_vector_3d: vec_xyz = vu.to_slope_plane(slope_vector_2d, slope_vector_2d)
    return ProjectedVectors(wind_vector_3d, slope_vector_3d)


@cy.cfunc
@cy.exceptval(check=False)
def get_phi_E(wind_vector_3d: vec_xyz, slope_vector_3d: vec_xyz, phi_W: cy.float, phi_S: cy.float) -> vec_xyz:
    # Convert wind and slope vectors to unit vectors on the slope-tangential plane
    w_S: vec_xyz = vu.as_unit_vector_3d(wind_vector_3d)  if phi_W > 0.0 else wind_vector_3d
    u_S: vec_xyz = vu.as_unit_vector_3d(slope_vector_3d) if phi_S > 0.0 else slope_vector_3d
    # Create the 3D slope-tangential phi_W, phi_S, and phi_E vectors
    phi_W_3d: vec_xyz = vu.scale_3d(phi_W, w_S)
    phi_S_3d: vec_xyz = vu.scale_3d(phi_S, u_S)
    phi_E_3d: vec_xyz = vu.add_3d(phi_W_3d, phi_S_3d)
    return phi_E_3d
# surface-fire-combine-wind-and-slope-vectors ends here
# [[file:../../org/pyretechnics.org::surface-fire-eccentricity][surface-fire-eccentricity]]
# TODO: Change model from str to enumerated type
@cy.cfunc
def surface_length_to_width_ratio(effective_wind_speed: cy.float, model: str = "behave") -> cy.float:
    """
    Calculate the length_to_width_ratio of the surface fire front given:
    - effective_wind_speed :: m/min (aligned with the slope-tangential plane)
    - model                :: "rothermel" or "behave" (Optional)
    """
    effective_wind_speed_mph = conv.m_min_to_mph(effective_wind_speed)
    # Select formula by model
    if model == "rothermel":
        return 1.0 + 0.25 * effective_wind_speed_mph

    elif model == "behave":
        return min(8.0,
                   0.936 * exp(0.1147 * effective_wind_speed_mph)
                   +
                   0.461 * exp(-0.0692 * effective_wind_speed_mph)
                   -
                   0.397)

    else:
        raise ValueError("Invalid input: model must be 'rothermel' or 'behave'.")


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def surface_fire_eccentricity(length_to_width_ratio: cy.float) -> cy.float:
    """
    Calculate the eccentricity (E) of the surface fire front using eq. 8 from
    Albini and Chase 1980 given:
    - L/W :: (1: circular spread, > 1: elliptical spread)
    """
    return sqrt((length_to_width_ratio * length_to_width_ratio) - 1.0) / length_to_width_ratio
# surface-fire-eccentricity ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-max][surface-fire-behavior-max]]
@cy.cfunc
@cy.exceptval(check=False)
def maybe_limit_wind_speed(use_wind_limit  : cy.bint,
                           max_wind_speed  : cy.float,
                           surface_fire_min: FireBehaviorMin,
                           phi_E_magnitude : cy.float) -> tuple[cy.float, cy.float]:
    """
    Given these inputs:
    - use_wind_limit   :: boolean
    - max_wind_speed   :: m/min
    - surface_fire_min :: FireBehaviorMin struct of no-wind-no-slope surface fire behavior values
      - base_spread_rate         :: m/min
      - base_fireline_intensity  :: kW/m
      - max_effective_wind_speed :: m/min
      - _phiS_G                  :: intermediate value for computing phi_S (unitless)
      - _phiW_scalr              :: intermediate value for computing phi_W (unitless)
      - _phiW_expnt              :: intermediate value for computing phi_W (unitless)
      - _ws_scalr                :: intermediate value for computing effective_wind_speed (m/min)
      - _ws_expnt                :: intermediate value for computing effective_wind_speed (unitless)
    - phi_E_magnitude  :: unitless

    return a tuple with these fields:
    - limited_wind_speed :: m/min
    - limited_phi_E      :: unitless
    """
    effective_wind_speed: cy.float = get_wind_speed(surface_fire_min, phi_E_magnitude)
    if (use_wind_limit and effective_wind_speed > max_wind_speed):
        return (
            max_wind_speed,
            get_phi_W(surface_fire_min, max_wind_speed),
        )
    else:
        return (
            effective_wind_speed,
            phi_E_magnitude,
        )


@cy.ccall
@cy.exceptval(check=False)
def calc_surface_fire_behavior_max(surface_fire_min      : FireBehaviorMin,
                                   midflame_wind_speed   : cy.float,
                                   upwind_direction      : cy.float,
                                   slope                 : cy.float,
                                   aspect                : cy.float,
                                   use_wind_limit        : cy.bint = True,
                                   surface_lw_ratio_model: str = "behave") -> FireBehaviorMax:
    """
    Given these inputs:
    - surface_fire_min       :: FireBehaviorMin struct of no-wind-no-slope surface fire behavior values
      - base_spread_rate         :: m/min
      - base_fireline_intensity  :: kW/m
      - max_effective_wind_speed :: m/min
      - _phiS_G                  :: intermediate value for computing phi_S (unitless)
      - _phiW_scalr              :: intermediate value for computing phi_W (unitless)
      - _phiW_expnt              :: intermediate value for computing phi_W (unitless)
      - _ws_scalr                :: intermediate value for computing effective_wind_speed (m/min)
      - _ws_expnt                :: intermediate value for computing effective_wind_speed (unitless)
    - midflame_wind_speed    :: m/min
    - upwind_direction       :: degrees clockwise from North
    - slope                  :: rise/run
    - aspect                 :: degrees clockwise from North
    - use_wind_limit         :: boolean (Optional)
    - surface_lw_ratio_model :: "rothermel" or "behave" (Optional)

    return a FireBehaviorMax struct containing these keys:
    - max_fire_type          :: 1 (surface)
    - max_spread_rate        :: m/min
    - max_spread_direction   :: (x, y, z) unit vector
    - max_fireline_intensity :: kW/m
    - max_flame_length       :: m
    - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
    - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
    - critical_spread_rate   :: m/min
    """
    # Unpack no-wind-no-slope surface fire behavior values
    spread_rate       : cy.float = surface_fire_min.base_spread_rate
    fireline_intensity: cy.float = surface_fire_min.base_fireline_intensity
    max_wind_speed    : cy.float = surface_fire_min.max_effective_wind_speed
    # Reverse the provided wind and slope directions
    downwind_direction: cy.float = conv.opposite_direction(upwind_direction)
    upslope_direction : cy.float = conv.opposite_direction(aspect)
    # Project wind and slope vectors onto the slope-tangential plane
    vectors        : ProjectedVectors = project_wind_and_slope_vectors_3d(midflame_wind_speed,
                                                                          downwind_direction,
                                                                          slope,
                                                                          upslope_direction)
    wind_vector_3d : vec_xyz = vectors.wind_vector_3d  # m/min
    slope_vector_3d: vec_xyz = vectors.slope_vector_3d # rise/run
    # Calculate phi_W and phi_S
    # NOTE: |wind_vector_3d| = slope-aligned midflame wind speed
    phi_W: cy.float = get_phi_W(surface_fire_min, vu.vector_magnitude_3d(wind_vector_3d))
    phi_S: cy.float = get_phi_S(surface_fire_min, slope)
    # Calculate phi_E and the max_spread_direction
    phi_E_3d: vec_xyz  = get_phi_E(wind_vector_3d, slope_vector_3d, phi_W, phi_S)
    phi_E   : cy.float = vu.vector_magnitude_3d(phi_E_3d)
    max_spread_direction: vec_xyz
    if phi_E > 0.0:
        max_spread_direction = vu.as_unit_vector_3d(phi_E_3d)
    elif phi_S > 0.0:
        max_spread_direction = vu.as_unit_vector_3d(slope_vector_3d)
    else:
        max_spread_direction = (0.0, 1.0, 0.0) # default: North
    # Limit effective wind speed to max wind speed if use_wind_limit == True
    (limited_wind_speed, limited_phi_E) = maybe_limit_wind_speed(use_wind_limit,
                                                                 max_wind_speed,
                                                                 surface_fire_min,
                                                                 phi_E)
    # Calculate and return max surface fire behavior values
    max_spread_rate       : cy.float = spread_rate * (1.0 + limited_phi_E)
    max_fireline_intensity: cy.float = fireline_intensity * (1.0 + limited_phi_E)
    length_to_width_ratio : cy.float = surface_length_to_width_ratio(limited_wind_speed, surface_lw_ratio_model)
    return FireBehaviorMax(
        max_fire_type          = 1,
        max_spread_rate        = max_spread_rate,
        max_spread_direction   = max_spread_direction,
        max_fireline_intensity = max_fireline_intensity,
        max_flame_length       = calc_flame_length(max_fireline_intensity),
        length_to_width_ratio  = length_to_width_ratio,
        eccentricity           = surface_fire_eccentricity(length_to_width_ratio),
        critical_spread_rate   = 0.0,
    )
# surface-fire-behavior-max ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-in-direction][surface-fire-behavior-in-direction]]
@cy.ccall
@cy.exceptval(check=False)
def calc_surface_fire_behavior_in_direction(surface_fire_max: FireBehaviorMax,
                                            spread_direction: vec_xyz) -> SpreadBehavior:
    """
    Given these inputs:
    - surface_fire_max     :: FireBehaviorMax struct of max surface fire behavior values
      - max_fire_type          :: 1 (surface)
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_fireline_intensity :: kW/m
      - max_flame_length       :: m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
      - critical_spread_rate   :: m/min
    - spread_direction     :: 3D unit vector on the slope-tangential plane

    return a SpreadBehavior struct containing these keys:
    - dphi_dt            :: phi/min
    - fire_type          :: 1 (surface)
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    # Unpack max surface fire behavior values
    max_spread_rate       : cy.float = surface_fire_max.max_spread_rate
    max_spread_direction  : vec_xyz  = surface_fire_max.max_spread_direction
    max_fireline_intensity: cy.float = surface_fire_max.max_fireline_intensity
    eccentricity          : cy.float = surface_fire_max.eccentricity
    # Calculate cos(w), where w is the offset angle between these unit vectors on the slope-tangential plane
    cos_w: cy.float = vu.dot_3d(max_spread_direction, spread_direction)
    # Calculate adjustment due to the offset angle from the max spread direction
    adjustment: cy.float = (1.0 - eccentricity) / (1.0 - eccentricity * cos_w)
    # Update surface fire behavior values by the adjustment value
    fireline_intensity: cy.float = max_fireline_intensity * adjustment
    return SpreadBehavior(
        dphi_dt            = 0.0,
        fire_type          = 1, # surface
        spread_rate        = max_spread_rate * adjustment,
        spread_direction   = spread_direction,
        fireline_intensity = fireline_intensity,
        flame_length       = calc_flame_length(fireline_intensity),
    )
# surface-fire-behavior-in-direction ends here
