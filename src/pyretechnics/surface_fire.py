# [[file:../../org/pyretechnics.org::surface-fire-imports][surface-fire-imports]]
# cython: profile=True
import cython
if cython.compiled:
    from cython.cimports.pyretechnics.math import exp, log, sqrt
    from cython.cimports.pyretechnics.cy_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, FireBehaviorMax
    import cython.cimports.pyretechnics.vector_utils as vu
else:
    from math import exp, log, sqrt
    from pyretechnics.py_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, FireBehaviorMax
    import pyretechnics.vector_utils as vu


import cython as cy

from pyretechnics.fuel_models import map_category, map_size_class, category_sum, size_class_sum
# surface-fire-imports ends here
# [[file:../../org/pyretechnics.org::surface-fire-common-intermediate-calculations][surface-fire-common-intermediate-calculations]]
def calc_surface_area_to_volume_ratio(f_i, f_ij, sigma):
    sigma_prime_i = size_class_sum(lambda i: f_ij[i] * sigma[i])
    return category_sum(lambda i: f_i[i] * sigma_prime_i[i])


def calc_packing_ratio(w_o, rho_p, delta):
    if (delta > 0.0):
        def scs_beta_i(i):
            return w_o[i] / rho_p[i]
        beta_i = size_class_sum(scs_beta_i) # TODO seems quite wasteful to realize this array
        def cs_beta_i(i):
            return beta_i[i]
        return category_sum(cs_beta_i) / delta
    else:
        return 0.0


def calc_optimum_packing_ratio(sigma_prime):
    return (3.348 / sigma_prime ** 0.8189) if (sigma_prime > 0.0) else 1.0
# surface-fire-common-intermediate-calculations ends here
# [[file:../../org/pyretechnics.org::surface-fire-reaction-intensity][surface-fire-reaction-intensity]]
def calc_mineral_damping_coefficients(f_ij, S_e):
    def scs_S_e_i(i):
        return f_ij[i] * S_e[i]
    S_e_i = size_class_sum(scs_S_e_i)
    def cs_mineral_damping_coefficient(i):
        S = S_e_i[i]
        return 0.174 / (S ** 0.19) if (S > 0.0) else 1.0
    return map_category(cs_mineral_damping_coefficient)


def calc_moisture_damping_coefficients(f_ij, M_f, M_x):
    def scs_M_f_i(i):
        return f_ij[i] * M_f[i]
    M_f_i = size_class_sum(scs_M_f_i)
    def scs_M_x_i(i):
        return f_ij[i] * M_x[i]
    M_x_i = size_class_sum(scs_M_x_i)
    def cs_moisture_damping_coefficient(i):
        M_f = M_f_i[i]
        M_x = M_x_i[i]
        r_M = min(1.0, M_f / M_x) if (M_x > 0.0) else 1.0
        return 1.0 - 2.59 * r_M + 5.11 * r_M ** 2.0 - 3.52 * r_M ** 3.0
    return map_category(cs_moisture_damping_coefficient)


def calc_low_heat_content(f_ij, h):
    def scs_h_i(i):
        return f_ij[i] * h[i]
    return size_class_sum(scs_h_i)


def calc_net_fuel_loading(g_ij, w_o, S_T):
    def scs_net_fuel_loading_i(i):
        g_ij_i = g_ij[i]
        w_o_i = w_o[i]
        S_T_i = S_T[i]
        return g_ij_i * w_o_i * (1.0 - S_T_i)
    return size_class_sum(scs_net_fuel_loading_i)


def calc_heat_per_unit_area(eta_S_i, eta_M_i, h_i, W_n_i):
    def cs_heat_per_unit_area(i):
        return W_n_i[i] * h_i[i] * eta_M_i[i] * eta_S_i[i]
    return category_sum(cs_heat_per_unit_area)


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
        def scs_rho_b_i(i):
            return w_o[i]
        rho_b_i = size_class_sum(scs_rho_b_i)
        def cs_ovendry_bulk_density(i):
            return rho_b_i[i]
        return category_sum(cs_ovendry_bulk_density) / delta
    else:
        return 0.0
# surface-fire-oven-dry-fuel-bed-bulk-density ends here
# [[file:../../org/pyretechnics.org::surface-fire-effective-heating-number-distribution][surface-fire-effective-heating-number-distribution]]
def calc_effective_heating_number_distribution(sigma):
    def msc_heating_number(i):
        sigma_i = sigma[i]
        return exp(-138.0 / sigma_i) if (sigma_i > 0.0) else 0.0
    return map_size_class(msc_heating_number)
# surface-fire-effective-heating-number-distribution ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-of-preignition-distribution][surface-fire-heat-of-preignition-distribution]]
def calc_heat_of_preignition_distribution(M_f):
    def msc_heat_of_preignition(i):
        M_f_i = M_f[i]
        return 250.0 + 1116.0 * M_f_i if (M_f_i > 0.0) else 0.0
    return map_size_class(msc_heat_of_preignition)
# surface-fire-heat-of-preignition-distribution ends here
# [[file:../../org/pyretechnics.org::surface-fire-heat-sink][surface-fire-heat-sink]]
def calc_heat_sink(f_i, f_ij, rho_b, epsilon_ij, Q_ig_ij):
    def scs_effective_heat_of_preignition_i(i):
        return f_ij[i] * epsilon_ij[i] * Q_ig_ij[i]
    effective_heat_of_preignition_i = size_class_sum(scs_effective_heat_of_preignition_i)
    def cs_effective_heat_of_preignition(i):
        return f_i[i] * effective_heat_of_preignition_i[i]
    effective_heat_of_preignition   = category_sum(cs_effective_heat_of_preignition)
    return rho_b * effective_heat_of_preignition
# surface-fire-heat-sink ends here
# [[file:../../org/pyretechnics.org::surface-fire-spread-rate-no-wind-no-slope][surface-fire-spread-rate-no-wind-no-slope]]
def calc_spread_rate(heat_source, heat_sink):
    return heat_source / heat_sink if (heat_sink > 0.0) else 0.0
# surface-fire-spread-rate-no-wind-no-slope ends here
# [[file:../../org/pyretechnics.org::surface-fire-intensity-functions][surface-fire-intensity-functions]]
def calc_residence_time(sigma_prime):
    """
    Returns the residence time (total burning time) of fuel (min) given:
    - sigma_prime :: ft^2/ft^3 (surface area to volume ratio)
    """
    return 384.0 / sigma_prime if (sigma_prime > 0.0) else 0.0


def calc_flame_depth(spread_rate, residence_time):
    """
    Returns the depth, or front-to-back distance, of the actively flaming zone
    of a free-spreading fire (ft) given:
    - spread_rate    :: ft/min (orthogonal to the fireline)
    - residence_time :: min
    """
    return spread_rate * residence_time


def calc_fireline_intensity(reaction_intensity, flame_depth):
    """
    Returns the rate of heat release per unit of fire edge (Btu/ft/s) given:
    - reaction_intensity :: Btu/ft^2/min
    - flame_depth        :: ft
    """
    return (reaction_intensity * flame_depth) / 60.0


def calc_flame_length(fireline_intensity):
    """
    Returns the average flame length (m) given:
    - fireline_intensity :: kW/m
    """
    return 0.07747042253266703 * (fireline_intensity ** 0.46)


def calc_areal_heat_output(spread_rate, fireline_intensity):
    """
    Returns the heat per unit area (kJ/m^2) given:
    - spread_rate        :: m/min
    - fireline_intensity :: kW/m
    """
    return 60.0 * fireline_intensity / spread_rate if spread_rate > 0.0 else 0.0
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


@cy.profile(False)
@cy.ccall
def calc_wind_adjustment_factor(fuel_bed_depth: cy.float, canopy_height: cy.float, canopy_cover: cy.float) -> cy.float:
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


#@cy.profile(False)
@cy.ccall
def calc_midflame_wind_speed(wind_speed_20ft: cy.float, fuel_bed_depth: cy.float, canopy_height: cy.float, canopy_cover: cy.float) -> cy.float:
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
import cython
if cython.compiled:
    from cython.cimports.pyretechnics.cy_types import vec_xy, vec_xyz
    from cython.cimports.pyretechnics.conversion import azimuthal_to_cartesian
    from cython.cimports.pyretechnics.vector_utils import as_unit_vector_3d, to_slope_plane, scale_3d, add_3d
else:
    from pyretechnics.py_types import vec_xy, vec_xyz
    from pyretechnics.conversion import azimuthal_to_cartesian
    from pyretechnics.vector_utils import as_unit_vector_3d, to_slope_plane, scale_3d, add_3d


import cython as cy


@cy.ccall
def project_wind_and_slope_vectors_3d(wind_speed: cy.float, downwind_direction: cy.float,
                                      slope: cy.float, upslope_direction: cy.float) -> dict[str,vec_xyz]:
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
    wind_vector_2d : vec_xy = azimuthal_to_cartesian(wind_speed, downwind_direction)
    slope_vector_2d: vec_xy = azimuthal_to_cartesian(slope, upslope_direction)
    # Project wind and slope vectors onto the slope-tangential plane
    wind_vector_3d : vec_xyz = to_slope_plane(wind_vector_2d, slope_vector_2d)
    slope_vector_3d: vec_xyz = to_slope_plane(slope_vector_2d, slope_vector_2d)
    return {
        "wind_vector_3d" : wind_vector_3d,
        "slope_vector_3d": slope_vector_3d,
    }


@cy.ccall
def get_phi_E(wind_vector_3d: vec_xyz, slope_vector_3d: vec_xyz, phi_W: cy.float, phi_S: cy.float) -> vec_xyz:
    # Convert wind and slope vectors to unit vectors on the slope-tangential plane
    w_S: vec_xyz = as_unit_vector_3d(wind_vector_3d)  if phi_W > 0.0 else wind_vector_3d
    u_S: vec_xyz = as_unit_vector_3d(slope_vector_3d) if phi_S > 0.0 else slope_vector_3d
    # Create the 3D slope-tangential phi_W, phi_S, and phi_E vectors
    phi_W_3d: vec_xyz = scale_3d(phi_W, w_S)
    phi_S_3d: vec_xyz = scale_3d(phi_S, u_S)
    phi_E_3d: vec_xyz = add_3d(phi_W_3d, phi_S_3d)
    return phi_E_3d
# surface-fire-combine-wind-and-slope-vectors ends here
# [[file:../../org/pyretechnics.org::surface-fire-eccentricity][surface-fire-eccentricity]]
from pyretechnics.conversion import m_min_to_mph


def surface_length_to_width_ratio(effective_wind_speed, model="rothermel"):
    """
    Calculate the length_to_width_ratio of the surface fire front given:
    - effective_wind_speed :: m/min (aligned with the slope-tangential plane)
    - model                :: "rothermel" or "behave" (Optional)
    """
    effective_wind_speed_mph = m_min_to_mph(effective_wind_speed)
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


def surface_fire_eccentricity(length_to_width_ratio):
    """
    Calculate the eccentricity (E) of the surface fire front using eq. 8 from
    Albini and Chase 1980 given:
    - L/W :: (1: circular spread, > 1: elliptical spread)
    """
    return sqrt(length_to_width_ratio ** 2.0 - 1.0) / length_to_width_ratio
# surface-fire-eccentricity ends here
# [[file:../../org/pyretechnics.org::surface-fire-behavior-max][surface-fire-behavior-max]]
import numpy as np
from pyretechnics.conversion import opposite_direction
from pyretechnics.vector_utils import vector_magnitude_3d


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
@cy.profile(True)
@cy.ccall
def calc_surface_fire_behavior_max(surface_fire_min: object, midflame_wind_speed: cy.float, upwind_direction: cy.float,
                                   slope: cy.float, aspect: cy.float,
                                   use_wind_limit: cy.bint,# = True, # FIXME optional can't seem to work in Cython, getting puzzling errors unexplained by documentation.
                                   surface_lw_ratio_model: object# = "rothermel"
                                   ) -> FireBehaviorMax:
    """
    Given these inputs:
    - surface_fire_min            :: dictionary of no-wind-no-slope surface fire behavior values
      - base_spread_rate             :: m/min
      - base_fireline_intensity      :: kW/m
      - max_effective_wind_speed     :: m/min
      - get_phi_S                    :: lambda: slope (rise/run) => phi_S (unitless)
      - get_phi_W                    :: lambda: midflame_wind_speed (m/min) => phi_W (unitless)
      - get_wind_speed               :: lambda: phi_W (unitless) => midflame_wind_speed (m/min)
    - midflame_wind_speed         :: m/min
    - upwind_direction            :: degrees clockwise from North
    - slope                       :: rise/run
    - aspect                      :: degrees clockwise from North
    - use_wind_limit              :: boolean (Optional)
    - surface_lw_ratio_model      :: "rothermel" or "behave" (Optional)

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
    phi_W = get_phi_W(vu.vector_magnitude_3d(wind_vector_3d)) # |wind_vector_3d| = slope-aligned midflame wind speed
    phi_S = get_phi_S(slope)
    # Calculate phi_E and the max_spread_direction
    phi_E_3d: vec_xyz = get_phi_E(wind_vector_3d, slope_vector_3d, phi_W, phi_S)
    phi_E: cy.float = vu.vector_magnitude_3d(phi_E_3d)
    max_spread_direction: vec_xyz
    if phi_E > 0.0:
        max_spread_direction = vu.as_unit_vector_3d(phi_E_3d)
    elif phi_S > 0.0:
        max_spread_direction = vu.as_unit_vector_3d(slope_vector_3d)
    else:
        max_spread_direction = (0.0, 1.0, 0.0) # default: North
    # Limit effective wind speed to max wind speed if use_wind_limit == True
    (limited_wind_speed, limited_phi_E) = maybe_limit_wind_speed(use_wind_limit, max_wind_speed,
                                                                 get_phi_W, get_wind_speed, phi_E)
    # Calculate and return max surface fire behavior values
    max_spread_rate        = spread_rate * (1.0 + limited_phi_E)
    max_fireline_intensity = fireline_intensity * (1.0 + limited_phi_E)
    length_to_width_ratio  = surface_length_to_width_ratio(limited_wind_speed, surface_lw_ratio_model)
    return FireBehaviorMax(
        -1.0, # FIXME
        max_spread_rate,
        max_spread_direction,
        max_fireline_intensity,
        calc_flame_length(max_fireline_intensity),
        length_to_width_ratio,
        surface_fire_eccentricity(length_to_width_ratio),
        -1.0, # FIXME
    )
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
    cos_w = np.dot(max_spread_direction, spread_direction)
    # Calculate adjustment due to the offset angle from the max spread direction
    adjustment = (1.0 - eccentricity) / (1.0 - eccentricity * cos_w)
    # Update surface fire behavior values by the adjustment value
    fireline_intensity = max_fireline_intensity * adjustment
    return {
        "fire_type"         : "surface",
        "spread_rate"       : max_spread_rate * adjustment,
        "spread_direction"  : np.asarray(spread_direction),
        "fireline_intensity": fireline_intensity,
        "flame_length"      : calc_flame_length(fireline_intensity),
    }
# surface-fire-behavior-in-direction ends here
