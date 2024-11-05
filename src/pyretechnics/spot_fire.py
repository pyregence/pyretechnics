# [[file:../../org/pyretechnics.org::expected-firebrand-production][expected-firebrand-production]]
from math import sqrt
import pyretechnics.surface_fire as sf


def expected_firebrand_production(fire_behavior, elevation_gradient, cube_resolution, firebrands_per_unit_heat=1e-6):
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
    if fire_behavior["spread_rate"] == 0.0:
        return 0.0
    else:
        #================================================================================================
        # Calculate the heat output per unit area
        #================================================================================================

        spread_rate          = fire_behavior["spread_rate"]                               # m/min
        fireline_intensity   = fire_behavior["fireline_intensity"]                        # kW/m
        heat_output_per_area = sf.calc_areal_heat_output(spread_rate, fireline_intensity) # kJ/m^2

        #================================================================================================
        # Calculate the slope-adjusted cell area
        #================================================================================================

        (dz_dx, dz_dy) = elevation_gradient                      # (rise/run, rise/run)
        slope_factor   = sqrt(1.0 + dz_dx ** 2.0 + dz_dy ** 2.0) # unitless
        cell_height    = cube_resolution[1]                      # meters
        cell_width     = cube_resolution[2]                      # meters
        cell_area      = cell_height * cell_width * slope_factor # m^2

        #================================================================================================
        # Calculate the expected number of firebrands produced in this cell
        #================================================================================================

        cell_heat_output = heat_output_per_area * cell_area            # kJ
        firebrand_count  = firebrands_per_unit_heat * cell_heat_output # number of firebrands
        return firebrand_count
# expected-firebrand-production ends here
# [[file:../../org/pyretechnics.org::convert-deltas][convert-deltas]]
def delta_to_grid_dx(cos_wdir, sin_wdir, delta_x, delta_y):
    """
    Computes the grid-aligned x coordinate of the delta vector, given the wind-aligned [ΔX ΔY] coordinates.
    Returns a signed distance (same unit as ΔX and ΔY).
    """
    wdir_x      = sin_wdir
    wdir_perp_x = cos_wdir
    return delta_x * wdir_x + delta_y * wdir_perp_x


def delta_to_grid_dy(cos_wdir, sin_wdir, delta_x, delta_y):
    """
    Computes the grid-aligned y coordinate of the delta vector, given the wind-aligned [ΔX ΔY] coordinates.
    Returns a signed distance (same unit as ΔX and ΔY).
    """
    wdir_y      = cos_wdir  # FIXME: Should this be negative or positive?
    wdir_perp_y = sin_wdir
    return delta_x * wdir_y + delta_y * wdir_perp_y


def distance_to_n_cells(distance, cell_size):
    """
    Converts a delta expressed as a signed distance to one expressed as a number of grid cells.
    """
    return round(distance / cell_size)
# convert-deltas ends here
# [[file:../../org/pyretechnics.org::resolve-spotting-lognormal-elmfire][resolve-spotting-lognormal-elmfire]]
from math import log, sqrt


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
from math import exp, sqrt, log
import pyretechnics.conversion as conv


def sample_normal(rand_gen, mu, sd):
    """
    Returns sample from normal/gaussian distribution given mu and sd.
    """
    return rand_gen.normal(mu, sd)


def sample_lognormal(rand_gen, mu, sd):
    """
    Returns sample from log-normal distribution given mu and sd.
    """
    return rand_gen.lognormal(mu, sd)


# FIXME: unused
def deltax_expected_value(mu_x, sigma_x):
    return conv.m_to_ft(exp(mu_x + (sigma_x ** 2.0) / 2.0))


# FIXME: unused
def deltax_coefficient_of_variation(sigma_x):
    return sqrt(exp(sigma_x ** 2.0) - 1.0)


def delta_x_sampler(spot_config, fireline_intensity, wind_speed_20ft):
    """
    Returns a function for randomly sampling ΔX, the spotting jump along the wind direction (in meters).
    """
    ln_params = resolve_lognormal_params(spot_config, fireline_intensity, wind_speed_20ft)
    mu_x      = ln_params["prob.lognormal.mu"]    # meters
    sigma_x   = ln_params["prob.lognormal.sigma"] # meters
    return lambda rand_gen: sample_lognormal(rand_gen, mu_x, sigma_x) # meters


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


def delta_y_sampler(spot_config, fireline_intensity, wind_speed_20ft):
    """
    Returns a function for randomly sampling ΔY, the spotting jump perpendicular to the wind direction (in meters).
    """
    sigma_y = resolve_crosswind_distance_stdev(spot_config, fireline_intensity, wind_speed_20ft) # meters
    return lambda rand_gen: sample_normal(rand_gen, 0.0, sigma_y) # meters
# sardoy-firebrand-dispersal ends here
# [[file:../../org/pyretechnics.org::sample-number-of-firebrands][sample-number-of-firebrands]]
def sample_poisson(rand_gen, mu):
    """
    Returns sample from poisson distribution given mu.
    """
    return rand_gen.poisson(mu)


def sample_number_of_firebrands(rand_gen, expected_firebrand_count):
    return sample_poisson(rand_gen, expected_firebrand_count)
# sample-number-of-firebrands ends here
# [[file:../../org/pyretechnics.org::firebrand-ignition-probability][firebrand-ignition-probability]]
from math import exp


def firebrand_flight_survival_probability(spotting_distance, decay_distance):
    """
    Returns the probability that a firebrand will survive its flight (Perryman 2012) given:
    - spotting_distance :: meters (d)
    - decay_distance    :: meters (1/lambda)

    P(Survival) = exp(-d * lambda)
    """
    return exp(-spotting_distance / decay_distance)


def heat_of_preignition(temperature, fine_fuel_moisture):
    """
    Returns heat of preignition given:
    - temperature        :: degrees Celsius
    - fine_fuel_moisture :: 0-1

    Q_ig = 144.512 - 0.266 * T_o - 0.00058 * (T_o)^2 - T_o * M + 18.54 * (1 - exp(-15.1 * M)) + 640 * M (eq. 10)
    """
    T_o = temperature
    M   = fine_fuel_moisture
    # Heat required to reach ignition temperature
    Q_a = 144.512 - 0.266 * T_o - 0.00058 * (T_o ** 2.0)
    # Heat required to raise moisture to the boiling point
    Q_b = -T_o * M
    # Heat of desorption
    Q_c = 18.54 * (1.0 - exp(-15.1 * M))
    # Heat required to vaporize moisture
    Q_d = 640.0 * M
    return Q_a + Q_b + Q_c + Q_d


def schroeder_ignition_probability(temperature, fine_fuel_moisture):
    """
    Returns the probability of spot fire ignition (Schroeder 1969) given:
    - temperature        :: degrees Celsius
    - fine_fuel_moisture :: 0-1

    X           = (400 - Q_ig) / 10
    P(Ignition) = (0.000048 * X^4.3) / 50 (pg. 15)
    """
    Q_ig        = heat_of_preignition(temperature, fine_fuel_moisture)
    X           = max(0.0, 400.0 - Q_ig) / 10.0
    P_Ignition  = 0.000048 * (X ** 4.3) / 50.0
    return min(P_Ignition, 1.0)
# firebrand-ignition-probability ends here
# [[file:../../org/pyretechnics.org::firebrands-time-of-ignition][firebrands-time-of-ignition]]
from math import sqrt
import pyretechnics.conversion as conv


def albini_firebrand_maximum_height(firebrand_diameter):
    return 0.39e5 * firebrand_diameter


def albini_t_max(flame_length):
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
    a           = 5.963  # dimensionless constant from (D33)
    b           = 4.563  # dimensionless constant from (D34)
    z           = 117.0  # maximum altitude of firebrands in meters [derived for (D44) in (Albini1979spot)]
    z_F         = flame_length                      # m
    w_F         = 2.3 * sqrt(flame_length)          # m/s
    charact_t   = conv.sec_to_min(2.0 * z_F / w_F)  # min
    # The following dimensionless factor is equal to t_T - t_o, with t_T defined by (D43) in Albini1979spot.
    travel_time = 1.2 + (a / 3.0) * (((b + z / z_F) / a) ** 1.5 - 1.0)
    return charact_t * travel_time


def spot_ignition_time(time_of_arrival, flame_length):
    """
    Returns the time of spot ignition using Albini 1979 and Perryman 2012 in minutes given:
    - time_of_arrival :: minutes
    - flame_length    :: meters

    t_spot = time_of_arrival + (2 * t_max) + t_ss
    """
    t_max          = albini_t_max(flame_length)
    t_steady_state = 20.0 # period of building up to steady state from ignition (min)
    return time_of_arrival + 2.0 * t_max + t_steady_state
# firebrands-time-of-ignition ends here
# [[file:../../org/pyretechnics.org::spread-firebrands][spread-firebrands]]
from math import sin, cos, hypot, radians
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.fuel_models as fm


def is_in_bounds(y, x, rows, cols):
    """
    Returns True if the grid coordinate (y,x) lies within the bounds [0,rows) by [0,cols).
    """
    return (y >= 0) and (x >= 0) and (y < rows) and (x < cols)


def is_burnable_cell(fuel_model_cube, t, y, x):
    """
    Returns True if the space-time coordinate (t,y,x) contains a burnable fuel model.
    """
    fuel_model_number = fuel_model_cube.get(t,y,x)
    fuel_model        = fm.fuel_model_table.get(fuel_model_number)
    return fuel_model and fuel_model["burnable"]


def cast_firebrand(rand_gen,
                   fuel_model_cube,
                   temperature_cube,
                   fuel_moisture_dead_1hr_cube,
                   fire_type_matrix,
                   firebrand_count_matrix, # NOTE: May be None
                   rows,
                   cols,
                   cell_height,
                   cell_width,
                   source_t,
                   source_y,
                   source_x,
                   decay_distance,
                   cos_wdir,
                   sin_wdir,
                   sample_delta_y_fn,
                   sample_delta_x_fn):
    """
    TODO: Add docstring
    Draws a random [ΔX, ΔY] pair of signed distances (in meters) from
    the supplied cell, representing the coordinates of the spotting jump in the directions
    parallel and perpendicular to the wind. ΔX will typically be positive (downwind),
    and positive ΔY means to the right of the downwind direction.
    """
    #=======================================================================================
    # Determine where the firebrand will land
    #=======================================================================================

    delta_y  = sample_delta_y_fn(rand_gen)                            # meters
    delta_x  = sample_delta_x_fn(rand_gen)                            # meters
    grid_dy  = delta_to_grid_dy(cos_wdir, sin_wdir, delta_x, delta_y) # meters
    grid_dx  = delta_to_grid_dx(cos_wdir, sin_wdir, delta_x, delta_y) # meters
    target_y = source_y + distance_to_n_cells(grid_dy, cell_height)
    target_x = source_x + distance_to_n_cells(grid_dx, cell_width)

    #=======================================================================================
    # Determine whether the firebrand will start a fire or fizzle out
    #=======================================================================================

    if is_in_bounds(target_y, target_x, rows, cols) and fire_type_matrix[target_y,target_x] == 0:
        # Firebrand landed on the grid in an unburned cell, so record it in firebrand_count_matrix (if provided)
        if isinstance(firebrand_count_matrix, np.ndarray):
            firebrand_count_matrix[target_y,target_x] += 1

        # Calculate the probability that the firebrand survived its flight and landed while still burning
        spotting_distance           = hypot(grid_dx, grid_dy) # meters
        flight_survival_probability = firebrand_flight_survival_probability(spotting_distance, decay_distance)

        # Roll the dice
        uniform_sample = rand_gen.uniform(0.0, 1.0)

        if (uniform_sample <= flight_survival_probability
            and is_burnable_cell(fuel_model_cube, source_t, target_y, target_x)):
            # Firebrand landed in a cell with a burnable fuel model, so calculate its ignition probability
            temperature          = temperature_cube.get(source_t, target_y, target_x)            # degrees Celsius
            fine_fuel_moisture   = fuel_moisture_dead_1hr_cube.get(source_t, target_y, target_x) # kg/kg
            ignition_probability = schroeder_ignition_probability(temperature, fine_fuel_moisture)

            if uniform_sample <= flight_survival_probability * ignition_probability:
                # Firebrand ignited the target cell, so return its coordinates for later processing
                return (target_y, target_x)


def spread_firebrands(space_time_cubes, output_matrices, cube_resolution, space_time_coordinate,
                      rand_gen, expected_firebrand_count, spot_config):
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
    - rand_gen                  :: numpy.random.Generator
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

    num_firebrands = sample_number_of_firebrands(rand_gen, expected_firebrand_count)

    if num_firebrands > 0:

        #=======================================================================================
        # Ensure that there is wind to transport the firebrands
        #=======================================================================================

        (t, y, x)      = space_time_coordinate
        wind_speed_10m = space_time_cubes["wind_speed_10m"].get(t,y,x) # km/hr

        if wind_speed_10m > 0.0:

            #=======================================================================================
            # Unpack all firebrand-related features of the source cell
            #=======================================================================================

            fuel_model_cube              = space_time_cubes["fuel_model"]
            temperature_cube             = space_time_cubes["temperature"]
            fuel_moisture_dead_1hr_cube  = space_time_cubes["fuel_moisture_dead_1hr"]
            fire_type_matrix             = output_matrices["fire_type"]
            firebrand_count_matrix       = output_matrices.get("firebrand_count")
            (_, rows, cols)              = fuel_model_cube.shape
            (_, cell_height, cell_width) = cube_resolution
            decay_distance               = spot_config["decay_distance"]
            upwind_direction             = space_time_cubes["upwind_direction"].get(t,y,x)
            downwind_direction           = radians(conv.opposite_direction(upwind_direction))
            cos_wdir                     = cos(downwind_direction)
            sin_wdir                     = sin(downwind_direction)
            fireline_intensity           = output_matrices["fireline_intensity"][y,x]             # m
            wind_speed_20ft              = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr
            wind_speed_20ft_mps          = conv.km_hr_to_mps(wind_speed_20ft)                     # m/s
            sample_delta_y_fn            = delta_y_sampler(spot_config, fireline_intensity, wind_speed_20ft_mps)
            sample_delta_x_fn            = delta_x_sampler(spot_config, fireline_intensity, wind_speed_20ft_mps)

            #=======================================================================================
            # Cast each firebrand, update firebrand_count_matrix, and accumulate any ignited cells
            #=======================================================================================

            ignited_cells = {ignited_cell for _i in range(num_firebrands)
                             if (ignited_cell := cast_firebrand(rand_gen,
                                                                fuel_model_cube,
                                                                temperature_cube,
                                                                fuel_moisture_dead_1hr_cube,
                                                                fire_type_matrix,
                                                                firebrand_count_matrix,
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
                                                                sample_delta_y_fn,
                                                                sample_delta_x_fn))
                             is not None}

            #=======================================================================================
            # Return any cells ignited by firebrands along with their time of ignition
            #=======================================================================================

            if len(ignited_cells) > 0:
                time_of_arrival = output_matrices["time_of_arrival"][y,x]           # minutes
                flame_length    = output_matrices["flame_length"][y,x]              # meters
                ignition_time   = spot_ignition_time(time_of_arrival, flame_length) # minutes

                return (ignition_time, ignited_cells)
# spread-firebrands ends here
