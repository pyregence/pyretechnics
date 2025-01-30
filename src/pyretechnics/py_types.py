# [[file:../../org/pyretechnics.org::py-types-py][py-types-py]]
#==============================================================
# Runtime-defined type aliases
#==============================================================

import cython as cy

pyidx     = cy.typedef(cy.Py_ssize_t)
vec_xy    = cy.typedef(tuple[cy.float, cy.float])
vec_xyz   = cy.typedef(tuple[cy.float, cy.float, cy.float])
coord_yx  = cy.typedef(tuple[pyidx, pyidx])
coord_tyx = cy.typedef(tuple[pyidx, pyidx, pyidx])
fcatarr   = cy.typedef(tuple[cy.float, cy.float])
fclaarr   = cy.typedef(tuple[cy.float, cy.float, cy.float, cy.float, cy.float, cy.float])

FuelModel = cy.struct(
    number               = cy.int,
    delta                = cy.float,
    M_x                  = fclaarr,
    M_f                  = fclaarr,
    w_o                  = fclaarr,
    sigma                = fclaarr,
    h                    = fclaarr,
    rho_p                = fclaarr,
    S_T                  = fclaarr,
    S_e                  = fclaarr,
    dynamic              = cy.bint,
    burnable             = cy.bint,
    exp_A_sigma          = fclaarr,
    firemod_size_classes = fclaarr,
    f_ij                 = fclaarr,
    f_i                  = fcatarr,
    g_ij                 = fclaarr,
)

ProjectedVectors = cy.struct(
    wind_vector_3d  = vec_xyz,
    slope_vector_3d = vec_xyz,
)

FireBehaviorMin = cy.struct(
    base_spread_rate         = cy.float,
    base_fireline_intensity  = cy.float,
    max_effective_wind_speed = cy.float,
    _phiS_G                  = cy.float,
    _phiW_scalr              = cy.float,
    _phiW_expnt              = cy.float,
    _ws_scalr                = cy.float,
    _ws_expnt                = cy.float,
)

FireBehaviorMax = cy.struct(
    max_fire_type          = cy.int,
    max_spread_rate        = cy.float,
    max_spread_direction   = vec_xyz,
    max_fireline_intensity = cy.float,
    max_flame_length       = cy.float,
    length_to_width_ratio  = cy.float,
    eccentricity           = cy.float,
    critical_spread_rate   = cy.float,
)

SpreadBehavior = cy.struct(
    dphi_dt            = cy.float,
    fire_type          = cy.int,
    spread_rate        = cy.float,
    spread_direction   = vec_xyz,
    fireline_intensity = cy.float,
    flame_length       = cy.float,
)

CrownSpreadInfo = cy.struct(
    fire_type            = cy.int,
    spread_rate          = cy.float,
    critical_spread_rate = cy.float,
)

SpotConfig = cy.struct(
    random_seed                  = cy.long,
    firebrands_per_unit_heat     = cy.float,
    downwind_distance_mean       = cy.float,
    fireline_intensity_exponent  = cy.float,
    wind_speed_exponent          = cy.float,
    downwind_variance_mean_ratio = cy.float,
    crosswind_distance_stdev     = cy.float,
    decay_distance               = cy.float,
)

# Pre-computed coefficients to apply elliptical wavelet math as fast as possible
# once the phi gradient information is available.
# See `pyretechnics.eulerian_level_set.dphi_dt_from_partialed_wavelet`.
PartialedEllWavelet = cy.struct(
    Vh_3d = vec_xyz,  # Heading spread rate vector (m/min)
    ewc_A = cy.float, # Dimensionless elliptical wavelet coefficient (<= 0)
    ewc_B = cy.float, # Dimensionless elliptical wavelet coefficient (<= 0)
    ewc_C = cy.float, # Dimensionless elliptical wavelet coefficient (>= 0)
)
# py-types-py ends here
