# [[file:../../org/pyretechnics.org::py-types][py-types]]
#==============================================================
# Runtime-defined type aliases
#==============================================================

import cython as cy

pyidx            = cy.typedef(cy.Py_ssize_t)
vec_xy           = cy.typedef(tuple) # tuple[cy.float, cy.float]
vec_xyz          = cy.typedef(tuple) # tuple[cy.float, cy.float, cy.float]
coord_yx         = cy.typedef(tuple) # tuple[pyidx, pyidx]
coord_tyx        = cy.typedef(tuple) # tuple[pyidx, pyidx, pyidx]
fcatarr          = cy.typedef(tuple) # tuple[cy.float, cy.float]
fclaarr          = cy.typedef(tuple) # tuple[cy.float, cy.float, cy.float, cy.float, cy.float, cy.float]
CompactFuelModel = cy.typedef(tuple) # tuple[cy.float, cy.float, cy.float, cy.float, cy.float, cy.float,
                                     #       cy.float, cy.float, cy.float, cy.float, cy.float, cy.float, cy.float]

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
    random_seed                  = cy.longlong,
    firebrands_per_unit_heat     = cy.float,
    downwind_distance_mean       = cy.float,
    fireline_intensity_exponent  = cy.float,
    wind_speed_exponent          = cy.float,
    downwind_variance_mean_ratio = cy.float,
    crosswind_distance_stdev     = cy.float,
    decay_distance               = cy.float,
)

JumpDistribution = cy.struct(
    # Downwind LogNormal params
    # Formally, we have ln(downwind_jump / 1m) ~ Normal(mu = mu_x, sigma = sigma_x)
    mu_x    = cy.float, # dimensionless (log-space)
    sigma_x = cy.float, # dimensionless (log-space)
    # Crosswind normal params
    # Formally, we have crosswind_jump ~ Normal(mu = 0, sigma = sigma_y)
    sigma_y = cy.float, # meters
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

CellInputs = cy.struct(
    slope                         = cy.float,
    aspect                        = cy.float,
    fuel_model_number             = cy.float,
    canopy_cover                  = cy.float,
    canopy_height                 = cy.float,
    canopy_base_height            = cy.float,
    canopy_bulk_density           = cy.float,
    wind_speed_10m                = cy.float,
    upwind_direction              = cy.float,
    fuel_moisture_dead_1hr        = cy.float,
    fuel_moisture_dead_10hr       = cy.float,
    fuel_moisture_dead_100hr      = cy.float,
    fuel_moisture_live_herbaceous = cy.float,
    fuel_moisture_live_woody      = cy.float,
    foliar_moisture               = cy.float,
    fuel_spread_adjustment        = cy.float,
    weather_spread_adjustment     = cy.float,
)

# Pre-computed information required to compute dphi/dt once the phi
# gradient is known. Derived from the surface and crown wavelets.
#
# NOTE: The reason to make this a small struct stored in an array is
#       efficiency - we want the CPU to have a low cache miss rate.
#
# NOTE: A significant benefit of this architecture is that it's
#       Rothermel-agnostic. EllipticalInfo could conceivably be
#       implemented using variants of the Rothermel model. This can be
#       valuable to give flexibility to users.
EllipticalInfo = cy.struct(
    cell_index           = coord_yx,
    elevation_gradient   = vec_xy,
    surface_wavelet      = PartialedEllWavelet,
    crown_wavelet        = PartialedEllWavelet,
    crowning_spread_rate = cy.float, # Surface spread rate at which crowning occurs
)

# Some data saved during the 1st Runge-Kutta pass.
Pass1CellOutput = cy.struct(
    cell_index      = coord_yx,
    phi_gradient_xy = vec_xy,
    dphi_dt_flim    = cy.float, # Flux-limited dphi/dt (phi/min, <= 0).
    phi_old         = cy.float,
)
# py-types ends here
