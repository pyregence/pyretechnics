# [[file:../../org/pyretechnics.org::expected-firebrand-production][expected-firebrand-production]]
# cython: profile=True
import cython
if cython.compiled:
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.cy_types as py_types
    from cython.cimports.pyretechnics.random import BufferedRandGen
    from cython.cimports.pyretechnics.space_time_cube import ISpaceTimeCube
    import cython.cimports.pyretechnics.surface_fire1 as sf
    from cython.cimports.pyretechnics.math import exp, sqrt, log, sin, cos, pow, round
else:
    import pyretechnics.conversion as conv
    import pyretechnics.py_types as py_types
    from pyretechnics.random import BufferedRandGen
    from pyretechnics.space_time_cube import ISpaceTimeCube
    import pyretechnics.surface_fire1 as sf
    from math import exp, sqrt, log, sin, cos, pow, round
import numpy as np
import cython as cy




@cy.ccall
@cy.exceptval(-1)
def expected_firebrand_production(
        fire_behavior: py_types.SpreadBehavior,
        elevation_gradient: py_types.vec_xy, 
        cell_horizontal_area_m2: cy.float,
        firebrands_per_unit_heat: cy.float=1e-6
        ) -> cy.float:
    """
    Return the expected number of firebrands produced by an entire cell when it burns given:
    - fire_behavior            :: dictionary of surface or crown fire behavior values
      - fire_type                 :: "unburned", "surface", "passive_crown", or "active_crown"
      - spread_rate               :: m/min
      - spread_direction          :: (x, y, z) unit vector on the slope-tangential plane
      - fireline_intensity        :: kW/m
      - flame_length              :: m
    - elevation_gradient       :: tuple with these fields
      - dz_dx                     :: rise/run
      - dz_dy                     :: rise/run
    - cube_resolution          :: tuple with these fields
      - band_duration             :: minutes
      - cell_height               :: meters
      - cell_width                :: meters
    - firebrands_per_unit_heat :: firebrands/kJ
    """
    if fire_behavior.spread_rate == 0.0:
        return 0.0
    else:
        #================================================================================================
        # Calculate the heat output per unit area
        #================================================================================================

        spread_rate: cy.float          = fire_behavior.spread_rate                               # m/min
        fireline_intensity: cy.float   = fire_behavior.fireline_intensity                        # kW/m
        heat_output_per_area: cy.float = sf.calc_areal_heat_output(spread_rate, fireline_intensity) # kJ/m^2

        #================================================================================================
        # Calculate the slope-adjusted cell area
        #================================================================================================

        (dz_dx, dz_dy) = elevation_gradient                      # (rise/run, rise/run)
        slope_factor: cy.float   = sqrt(1.0 + (dz_dx * dz_dx) + (dz_dy * dz_dy)) # unitless
        cell_area: cy.float      = cell_horizontal_area_m2 * slope_factor # m^2

        #================================================================================================
        # Calculate the expected number of firebrands produced in this cell
        #================================================================================================

        cell_heat_output: cy.float = heat_output_per_area * cell_area            # kJ
        firebrand_count: cy.float  = firebrands_per_unit_heat * cell_heat_output # number of firebrands
        return firebrand_count
# expected-firebrand-production ends here
# [[file:../../org/pyretechnics.org::convert-deltas][convert-deltas]]
@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def delta_to_grid_dx(cos_wdir: cy.float, sin_wdir: cy.float, delta_x: cy.float, delta_y: cy.float) -> cy.float:
    """
    Computes the grid-aligned x coordinate of the delta vector, given the wind-aligned [ΔX ΔY] coordinates.
    Returns a signed distance (same unit as ΔX and ΔY).
    """
    wdir_x: cy.float      = sin_wdir
    wdir_perp_x: cy.float = cos_wdir
    return delta_x * wdir_x + delta_y * wdir_perp_x


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def delta_to_grid_dy(cos_wdir: cy.float, sin_wdir: cy.float, delta_x: cy.float, delta_y: cy.float) -> cy.float:
    """
    Computes the grid-aligned y coordinate of the delta vector, given the wind-aligned [ΔX ΔY] coordinates.
    Returns a signed distance (same unit as ΔX and ΔY).
    """
    wdir_y: cy.float      = cos_wdir  # FIXME: Should this be negative or positive?
    wdir_perp_y: cy.float = sin_wdir
    return delta_x * wdir_y + delta_y * wdir_perp_y


@cy.ccall
@cy.inline
@cy.cdivision(True)
def distance_to_n_cells(distance: cy.float, cell_size: cy.float) -> cy.float:
    """
    Converts a delta expressed as a signed distance to one expressed as a number of grid cells.
    """
    return round(distance / cell_size)
# convert-deltas ends here
# [[file:../../org/pyretechnics.org::resolve-spotting-lognormal-elmfire][resolve-spotting-lognormal-elmfire]]


def resolve_exp_delta_x(spot_config, fireline_intensity, wind_speed_20ft):
    """
    Computes the expected value E[ΔX] (in meters) of the downwind spotting distance ΔX given:
    - spot_config        :: a map of spotting parameters
    - fireline_intensity :: kW/m
    - wind_speed_20ft    :: m/s
    """
    downwind_distance_mean      = spot_config["downwind_distance_mean"]
    fireline_intensity_exponent = spot_config["fireline_intensity_exponent"]
    wind_speed_exponent         = spot_config["wind_speed_exponent"]
    return (downwind_distance_mean
            * (fireline_intensity ** fireline_intensity_exponent)
            * (wind_speed_20ft ** wind_speed_exponent))


def resolve_var_delta_x(spot_config, exp_delta_x):
    """
    Computes the variance Var[ΔX] (in m^2) of the downwind spotting distance ΔX given:
    - spot_config :: a map of spotting parameters
    - exp_delta_x :: meters (E[ΔX])
    """
    return spot_config["downwind_variance_mean_ratio"] * exp_delta_x


def lognormal_mu_from_moments(mean, variance):
    """
    TODO: Add docstring
    """
    m2 = mean ** 2.0
    return log(m2 / sqrt(m2 + variance))


def lognormal_sigma_from_moments(mean, variance):
    """
    TODO: Add docstring
    """
    return sqrt(log(1.0 + variance / (mean ** 2.0)))


def resolve_lognormal_params(spot_config, fireline_intensity, wind_speed_20ft):
    """
    TODO: Add docstring
    """
    exp_delta_x = resolve_exp_delta_x(spot_config, fireline_intensity, wind_speed_20ft)
    var_delta_x = resolve_var_delta_x(spot_config, exp_delta_x)
    return {
        "prob.lognormal.mu"   : lognormal_mu_from_moments(exp_delta_x, var_delta_x),
        "prob.lognormal.sigma": lognormal_sigma_from_moments(exp_delta_x, var_delta_x),
    }
# resolve-spotting-lognormal-elmfire ends here
# [[file:../../org/pyretechnics.org::sardoy-firebrand-dispersal][sardoy-firebrand-dispersal]]

@cy.ccall
def sample_normal(rng: BufferedRandGen, mu: cy.float, sd: cy.float) -> cy.float:
    """
    Returns sample from normal/gaussian distribution given mu and sd.
    """
    return mu + sd * (rng.next_normal())


@cy.ccall
def sample_lognormal(rng: BufferedRandGen, mu: cy.float, sd: cy.float) -> cy.float:
    """
    Returns sample from log-normal distribution given mu and sd.
    """
    return exp(sample_normal(rng, mu, sd))


# FIXME: unused
def deltax_expected_value(mu_x, sigma_x):
    return conv.m_to_ft(exp(mu_x + (sigma_x ** 2.0) / 2.0))


# FIXME: unused
def deltax_coefficient_of_variation(sigma_x):
    return sqrt(exp(sigma_x ** 2.0) - 1.0)


# When will we have the default sigma_Y > E[ΔX]?
# It can be seen that this nonsensical situation
# happens iff sigma_X exceeds the following number:
#
# sqrt(log(1.0 + (0.88 ** 2.0) / (0.92 * 0.47))
#
# => 1.0131023746492023
sigma_y_scalar_m = 0.92 * 0.47 / (0.88 ** 2.0)


def himoto_resolve_default_sigma_y_from_lognormal_params(mu_x, sigma_x):
    es2h       = exp((sigma_x ** 2.0) / 2.0)
    avg_deltax = exp(mu_x) * es2h
    return sigma_y_scalar_m * avg_deltax * (es2h + 1.0) * (es2h - 1.0) # meters


def himoto_resolve_default_sigma_y(spot_config, fireline_intensity, wind_speed_20ft):
    ln_params = resolve_lognormal_params(spot_config, fireline_intensity, wind_speed_20ft)
    mu_x      = ln_params["prob.lognormal.mu"]    # meters
    sigma_x   = ln_params["prob.lognormal.sigma"] # meters
    return himoto_resolve_default_sigma_y_from_lognormal_params(mu_x, sigma_x) # meters


def resolve_crosswind_distance_stdev(spot_config, fireline_intensity, wind_speed_20ft):
    crosswind_distance_stdev = spot_config.get("crosswind_distance_stdev")
    if crosswind_distance_stdev != None:
        return crosswind_distance_stdev # meters
    else:
        return himoto_resolve_default_sigma_y(spot_config, fireline_intensity, wind_speed_20ft) # meters


@cy.cclass
class JumpDistribution:
    # Downwind LogNormal params
    mu_x: cy.float # dimensionless (log-space)
    sigma_x: cy.float # dimensionless (log-space)
    # Formally, we have ln(downwind_jump / 1m) ~ Normal(mu = mu_x, sigma = sigma_x)
    
    # Cross-wind normal params
    sigma_y: cy.float # meters
    # Formally, we have crosswind_jump ~ Normal(mu = 0, sigma = sigma_y)

@cy.ccall
def resolve_JumpDistribution(spot_config, fireline_intensity: cy.float, wind_speed_20ft: cy.float) -> JumpDistribution:
    ret: JumpDistribution = JumpDistribution()
    ln_params = resolve_lognormal_params(spot_config, fireline_intensity, wind_speed_20ft)
    # FIXME why are these described as meters, when they're log-space quantities?
    ret.mu_x      = ln_params["prob.lognormal.mu"]    # meters
    ret.sigma_x   = ln_params["prob.lognormal.sigma"] # meters
    ret.sigma_y = resolve_crosswind_distance_stdev(spot_config, fireline_intensity, wind_speed_20ft) # meters
    return ret

@cy.ccall
def sample_downwind_jump(jd: JumpDistribution, random_generator: BufferedRandGen) -> cy.float:
    return sample_lognormal(random_generator, jd.mu_x, jd.sigma_x)

@cy.ccall
def sample_crosswind_jump(jd: JumpDistribution, random_generator: BufferedRandGen) -> cy.float:
    return sample_normal(random_generator, 0, jd.sigma_y)


# sardoy-firebrand-dispersal ends here
# [[file:../../org/pyretechnics.org::sample-number-of-firebrands][sample-number-of-firebrands]]
# sample-number-of-firebrands ends here
# [[file:../../org/pyretechnics.org::firebrand-ignition-probability][firebrand-ignition-probability]]


@cy.ccall
@cy.cdivision(True)
def firebrand_flight_survival_probability(spotting_distance: cy.float, decay_distance: cy.float) -> cy.float:
    """
    Returns the probability that a firebrand will survive its flight (Perryman 2012) given:
    - spotting_distance :: meters (d)
    - decay_distance    :: meters (1/lambda)

    P(Survival) = exp(-d * lambda)
    """
    return exp(-spotting_distance / decay_distance)

@cy.ccall
@cy.exceptval(check=False)
def heat_of_preignition(temperature: cy.float, fine_fuel_moisture: cy.float):
    """
    Returns heat of preignition given:
    - temperature        :: degrees Celsius
    - fine_fuel_moisture :: 0-1

    Q_ig = 144.512 - 0.266 * T_o - 0.00058 * (T_o)^2 - T_o * M + 18.54 * (1 - exp(-15.1 * M)) + 640 * M (eq. 10)
    """
    T_o: cy.float = temperature
    M: cy.float   = fine_fuel_moisture
    # Heat required to reach ignition temperature
    Q_a: cy.float = 144.512 - 0.266 * T_o - 0.00058 * (T_o*T_o)
    # Heat required to raise moisture to the boiling point
    Q_b: cy.float = -T_o * M
    # Heat of desorption
    Q_c: cy.float = 18.54 * (1.0 - exp(-15.1 * M))
    # Heat required to vaporize moisture
    Q_d: cy.float = 640.0 * M
    return Q_a + Q_b + Q_c + Q_d

@cy.ccall
@cy.exceptval(check=False)
def schroeder_ignition_probability(temperature: cy.float, fine_fuel_moisture: cy.float) -> cy.float:
    """
    Returns the probability of spot fire ignition (Schroeder 1969) given:
    - temperature        :: degrees Celsius
    - fine_fuel_moisture :: 0-1

    X           = (400 - Q_ig) / 10
    P(Ignition) = (0.000048 * X^4.3) / 50 (pg. 15)
    """
    Q_ig: cy.float        = heat_of_preignition(temperature, fine_fuel_moisture)
    X: cy.float           = max(0.0, 400.0 - Q_ig) * 0.1
    P_Ignition: cy.float  = 0.000048 * pow(X, 4.3) * 0.02
    return min(P_Ignition, 1.0)
# firebrand-ignition-probability ends here
# [[file:../../org/pyretechnics.org::firebrands-time-of-ignition][firebrands-time-of-ignition]]

@cy.ccall
@cy.exceptval(check=False)
@cy.inline
def albini_firebrand_maximum_height(firebrand_diameter: cy.float) -> cy.float:
    return 0.39e5 * firebrand_diameter

@cy.ccall
@cy.exceptval(-1)
@cy.cdivision(True)
def albini_t_max(flame_length: cy.float) -> cy.float:
    """
    Returns the time of spot ignition using Albini1979spot in minutes given:
    - flame_length :: meters [z_F]

    a           = 5.963                                                             (D33)
    b           = a - 1.4                                                           (D34)
    D           = 0.003
    z           = 0.39 * D * 10^5
    w_F         = 2.3 * z_F^0.5                                                     (A58)
    t_c         = 1
    t_o         = t_c / (2 * z_F / w_F)
    travel_time = t_1 + t_2 + t_3 = 1.2 + (a / 3) * (((b + (z/z_F)) / a)^3/2 - 1)   (D43)
    """
    a: cy.float           = 5.963  # dimensionless constant from (D33)
    b: cy.float           = 4.563  # dimensionless constant from (D34)
    z: cy.float           = 117.0  # maximum altitude of firebrands in meters [derived for (D44) in (Albini1979spot)]
    z_F: cy.float         = flame_length                      # m
    w_F: cy.float         = 2.3 * sqrt(flame_length)          # m/s
    charact_t: cy.float   = conv.sec_to_min(2.0 * z_F / w_F)  # min
    # The following dimensionless factor is equal to t_T - t_o, with t_T defined by (D43) in Albini1979spot.
    u: cy.float = (b + z / z_F) / a
    u3_2: cy.float = u*sqrt(u) # Faster than ** 1.5
    travel_time: cy.float = 1.2 + (a / 3.0) * (u3_2 - 1.0)
    return charact_t * travel_time

@cy.ccall
@cy.exceptval(-1)
@cy.cdivision(True)
def spot_ignition_time(time_of_arrival: cy.float, flame_length: cy.float) -> cy.float:
    """
    Returns the time of spot ignition using Albini 1979 and Perryman 2012 in minutes given:
    - time_of_arrival :: minutes
    - flame_length    :: meters

    t_spot = time_of_arrival + (2 * t_max) + t_ss
    """
    t_max: cy.float          = albini_t_max(flame_length) # FIXME what is this still doing here? I thought we decided to neglect the flight time! (Val)
    t_steady_state: cy.float = 20.0 # period of building up to steady state from ignition (min)
    return time_of_arrival + 2.0 * t_max + t_steady_state
# firebrands-time-of-ignition ends here
# [[file:../../org/pyretechnics.org::spread-firebrands][spread-firebrands]]
import numpy as np
import pyretechnics.fuel_models as fm


@cy.ccall
@cy.inline
def is_in_bounds(y: cy.int, x: cy.int, rows: py_types.pyidx, cols: py_types.pyidx) -> cy.bint:
    """
    Returns True if the grid coordinate (y,x) lies within the bounds [0,rows) by [0,cols).
    """
    return (y >= 0) and (x >= 0) and (y < rows) and (x < cols)

@cy.ccall
def is_burnable_cell(fuel_model_cube: ISpaceTimeCube, t: py_types.pyidx, y: py_types.pyidx, x: py_types.pyidx) -> cy.bint:
    """
    Returns True if the space-time coordinate (t,y,x) contains a burnable fuel model.
    """
    fuel_model_number: py_types.pyidx = cy.cast(py_types.pyidx, fuel_model_cube.get_float(t,y,x))
    fm_table: dict = fm.fuel_model_table
    fuel_model     = fm_table.get(fuel_model_number)  # FIXME INSANELY slow
    return fuel_model and fuel_model["burnable"]


@cy.ccall
def cast_firebrand(rng: BufferedRandGen,
                   fuel_model_cube: ISpaceTimeCube,
                   temperature_cube: ISpaceTimeCube,
                   fuel_moisture_dead_1hr_cube: ISpaceTimeCube,
                   firebrand_count_matrix, # NOTE: May be None
                   rows: py_types.pyidx,
                   cols: py_types.pyidx,
                   cell_height: cy.float,
                   cell_width: cy.float,
                   source_t: py_types.pyidx,
                   source_y: py_types.pyidx,
                   source_x: py_types.pyidx,
                   decay_distance: cy.float,
                   cos_wdir: cy.float,
                   sin_wdir: cy.float,
                   jd: JumpDistribution
                   ) -> py_types.coord_yx:
    """
    TODO: Add docstring
    Draws a random [ΔX, ΔY] pair of signed distances (in meters) from
    the supplied cell, representing the coordinates of the spotting jump in the directions
    parallel and perpendicular to the wind. ΔX will typically be positive (downwind),
    and positive ΔY means to the right of the downwind direction.

    Note that if the random draw yields no spot ignition, the source cell coordinates will be returned:
    calling code should check that.
    """
    #=======================================================================================
    # Determine where the firebrand will land
    #=======================================================================================

    delta_y: cy.float  = sample_crosswind_jump(jd, rng)                    # meters
    delta_x: cy.float  = sample_downwind_jump(jd, rng)                     # meters
    grid_dy: cy.float  = delta_to_grid_dy(cos_wdir, sin_wdir, delta_x, delta_y) # meters
    grid_dx: cy.float  = delta_to_grid_dx(cos_wdir, sin_wdir, delta_x, delta_y) # meters
    # NOTE it would cause a bug to type the following as py_types.pyidx.
    target_y: cy.int = source_y + cy.cast(cy.int, distance_to_n_cells(grid_dy, cell_height))
    target_x: cy.int = source_x + cy.cast(cy.int, distance_to_n_cells(grid_dx, cell_width))

    #=======================================================================================
    # Determine whether the firebrand will start a fire or fizzle out
    #=======================================================================================

    if is_in_bounds(target_y, target_x, rows, cols): # and fire_type_matrix[target_y,target_x] == 0: # NOTE REVIEW I removed this call because it is inconsistent. The target cell may well burn between firebrand emission and spot ignition.
        # Firebrand landed on the grid in an unburned cell, so record it in firebrand_count_matrix (if provided)
        if firebrand_count_matrix is not None:
            firebrand_count_matrix[target_y,target_x] += 1 # FIXME remove this useless performance hog.

        # Calculate the probability that the firebrand survived its flight and landed while still burning
        spotting_distance: cy.float           = sqrt(grid_dx*grid_dx + grid_dy*grid_dy) # meters
        flight_survival_probability: cy.float = firebrand_flight_survival_probability(spotting_distance, decay_distance)

        # Roll the dice
        uniform_sample: cy.float = rng.next_uniform()

        if (uniform_sample <= flight_survival_probability
            and is_burnable_cell(fuel_model_cube, source_t, target_y, target_x)):
            # Firebrand landed in a cell with a burnable fuel model, so calculate its ignition probability``
            temperature: cy.float          = temperature_cube.get_float(source_t, target_y, target_x)            # degrees Celsius
            fine_fuel_moisture: cy.float   = fuel_moisture_dead_1hr_cube.get_float(source_t, target_y, target_x) # kg/kg
            ignition_probability: cy.float = schroeder_ignition_probability(temperature, fine_fuel_moisture)

            if uniform_sample <= flight_survival_probability * ignition_probability:
                # Firebrand ignited the target cell, so return its coordinates for later processing
                return (target_y, target_x)
    # This code is only reached if the spotting ignition fails.
    return (source_y, source_x) # For efficiency, the source cell is used as a sentinel value for a failed spot ignition.

@cy.ccall
def spread_firebrands(
        fuel_model_cube: ISpaceTimeCube,
        temperature_cube: ISpaceTimeCube,
        fuel_moisture_dead_1hr_cube: ISpaceTimeCube,
        sim_area_bounds: py_types.coord_yx,
        spatial_resolution: py_types.vec_xy,
        space_time_coordinate: py_types.coord_tyx,
        upwind_direction: cy.float, # degrees
        wind_speed_10m: cy.float, # km/hr (for shame!)
        fireline_intensity: cy.float,
        flame_length: cy.float,
        time_of_arrival: cy.float,
        random_generator,
        num_firebrands: py_types.pyidx,
        spot_config):
    """
    Given these inputs:
    - space_time_cubes          :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - fuel_model                    :: integer index in fm.fuel_model_table
      - temperature                   :: degrees Celsius
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
    - output_matrices           :: dictionary of 2D Numpy arrays whose spatial dimensions match the space_time_cubes
      - fire_type                     :: 2D byte array (0=unburned, 1=surface, 2=passive_crown, 3=active_crown)
      - fireline_intensity            :: 2D float array (kW/m)
      - flame_length                  :: 2D float array (m)
      - time_of_arrival               :: 2D float array (min)
      - firebrand_count               :: 2D integer array (number of firebrands) (Optional)
    - cube_resolution           :: tuple with these fields
      - band_duration                 :: minutes
      - cell_height                   :: meters
      - cell_width                    :: meters
    - space_time_coordinate     :: (t,y,x) coordinate in which the source cell burns
    - random_generator          :: numpy.random.Generator
    - expected_firebrand_count  :: expected number of firebrands to cast
    - spot_config               :: dictionary of spotting parameters
      - downwind_distance_mean        :: meters
      - fireline_intensity_exponent   :: downwind_distance_mean multiplier [I^fireline_intensity_exponent]
      - wind_speed_exponent           :: downwind_distance_mean multiplier [U^wind_speed_exponent]
      - downwind_variance_mean_ratio  :: meters^2 / meter [downwind_variance_mean_ratio = Var(X) / E(X)]
      - crosswind_distance_stdev      :: meters
      - decay_distance                :: meters

    samples a number of firebrands from a Poisson distribution parameterized by expected_firebrand_count,
    casts these from the space_time_coordinate into grid cells in the space-time cube, records their landing
    locations in output_matrices["firebrand_count"] (if provided), filters out all of the firebrands that
    fizzle out in either burnable or non-burnable fuels, and returns any that ignite new spot fires in
    a tuple with these fields:

    - ignition_time :: minutes
    - ignited_cells :: set of (y,x) grid coordinates
    """
    #=======================================================================================
    # Sample the number of firebrands to cast from the source cell
    #=======================================================================================

    if num_firebrands > 0:

        #=======================================================================================
        # Ensure that there is wind to transport the firebrands
        #=======================================================================================

        (t, y, x)      = space_time_coordinate

        if wind_speed_10m > 0.0:

            #=======================================================================================
            # Unpack all firebrand-related features of the source cell
            #=======================================================================================
            firebrand_count_matrix       = None # FIXME
            (rows, cols)                 = sim_area_bounds
            (cell_height, cell_width)    = spatial_resolution
            decay_distance: cy.float               = spot_config["decay_distance"] # FIXME faster lookup (Extension Type)
            downwind_direction: cy.float           = conv.deg_to_rad(conv.opposite_direction(upwind_direction))
            # OPTIM get rid of trigonometry here if possible by having callers pass vectors.
            cos_wdir: cy.float                     = cos(downwind_direction)
            sin_wdir: cy.float                     = sin(downwind_direction)
            wind_speed_20ft: cy.float              = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr
            wind_speed_20ft_mps: cy.float          = conv.km_hr_to_mps(wind_speed_20ft)                     # m/s
            jd: JumpDistribution                   = resolve_JumpDistribution(spot_config, fireline_intensity, wind_speed_20ft_mps)

            #=======================================================================================
            # Cast each firebrand, update firebrand_count_matrix, and accumulate any ignited cells
            #=======================================================================================

            # FIXME a set is slow, use a list instead. Collisions can happen with other source cells anyway.
            ignited_cells = {ignited_cell for _i in range(num_firebrands) # FIXME native type.
                             if (ignited_cell := cast_firebrand(random_generator,
                                                                fuel_model_cube,
                                                                temperature_cube,
                                                                fuel_moisture_dead_1hr_cube,
                                                                firebrand_count_matrix, # FIXME clear
                                                                rows,
                                                                cols,
                                                                cell_height,
                                                                cell_width,
                                                                t,
                                                                y,
                                                                x,
                                                                decay_distance,
                                                                cos_wdir,
                                                                sin_wdir,
                                                                jd))
                             != (y, x)}

            #=======================================================================================
            # Return any cells ignited by firebrands along with their time of ignition
            #=======================================================================================
            if len(ignited_cells) > 0:
                ignition_time   = spot_ignition_time(time_of_arrival, flame_length) # minutes
                return (ignition_time, ignited_cells)
# spread-firebrands ends here
