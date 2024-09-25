# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients-approx][phi-field-spatial-gradients-approx]]
import numpy as np


def calc_dphi_dx_approx(phi, dx, x, y):
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given the cell width dx.

    NOTE: The origin cell (x=0,y=0) is located in the upper left corner
          of the grid in Python arrays. Thus, as x increases, we move
          to the east, and as y increases, we move to the south.
    """
    return (phi[y][x+1] - phi[y][x-1]) / (2.0 * dx)


def calc_dphi_dy_approx(phi, dy, x, y):
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given the cell height dy.

    NOTE: The origin cell (x=0,y=0) is located in the upper left corner
          of the grid in Python arrays. Thus, as x increases, we move
          to the east, and as y increases, we move to the south.
    """
    return (phi[y-1][x] - phi[y+1][x]) / (2.0 * dy)


def calc_phi_gradient_approx(phi, dx, dy, x, y):
    """
    Calculate the spatial gradient of the phi raster at grid cell (x,y)
    given the cell width dx and the cell height dy.
    """
    dphi_dx = calc_dphi_dx_approx(phi, dx, x, y)
    dphi_dy = calc_dphi_dy_approx(phi, dy, x, y)
    return np.asarray((dphi_dx, dphi_dy))
# phi-field-spatial-gradients-approx ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector][phi-field-normal-vector]]
import pyretechnics.vector_utils as vu


def calc_phi_normal_vector(phi, dx, dy, x, y):
    """
    Calculate the phi field normal vector in the x and y dimensions.

    - n_x: eastward component of the unit normal vector
    - n_y: northward component of the unit normal vector
    """
    phi_gradient = calc_phi_gradient_approx(phi, dx, dy, x, y)
    if phi_gradient[0] == 0.0 and phi_gradient[1] == 0.0:
        return phi_gradient # (n_x, n_y)
    else:
        return vu.as_unit_vector(phi_gradient) # (n_x, n_y)
# phi-field-normal-vector ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector-angle][phi-field-normal-vector-angle]]
from math import atan, pi, degrees


def calc_phi_normal_azimuth(phi_normal_vector):
    """
    Calculate the angle (measured in degrees clockwise from North)
    to which the phi field's normal vector points.
    """
    (n_x, n_y) = phi_normal_vector
    if n_x > 0:
        if n_y >= 0:
            angle = 1/2 * pi - atan(n_y / n_x)
        elif n_y < 0:
            angle = 1/2 * pi + atan(abs(n_y) / n_x)
    elif n_x < 0:
        if n_y >= 0:
            angle = 3/2 * pi + atan(n_y / abs(n_x))
        elif n_y < 0:
            angle = 3/2 * pi - atan(n_y / n_x)
    else:
        if n_y >= 0:
            angle = 0.0
        elif n_y < 0:
            angle = pi
    return degrees(angle)
# phi-field-normal-vector-angle ends here
# [[file:../../org/pyretechnics.org::superbee-flux-limiter][superbee-flux-limiter]]
def calc_superbee_flux_limiter(dphi_up, dphi_loc):
    """
    TODO: Add docstring
    """
    r = dphi_up / dphi_loc
    return max(0,
               min(2 * r, 1),
               min(r, 2))
# superbee-flux-limiter ends here
# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients][phi-field-spatial-gradients]]
import numpy as np


def calc_dphi_dx(phi, u_x, dx, x, y):
    """
    TODO: Add docstring
    """
    phi_east = calc_phi_east(phi, u_x, x, y)
    phi_west = calc_phi_west(phi, u_x, x, y)
    return (phi_east - phi_west) / dx


def calc_dphi_dy(phi, u_y, dy, x, y):
    """
    TODO: Add docstring
    """
    phi_north = calc_phi_north(phi, u_y, x, y)
    phi_south = calc_phi_south(phi, u_y, x, y)
    return (phi_north - phi_south) / dy


def calc_phi_gradient(phi, u_x, u_y, dx, dy, x, y):
    """
    TODO: Add docstring
    """
    dphi_dx = calc_dphi_dx(phi, u_x, dx, x, y)
    dphi_dy = calc_dphi_dy(phi, u_y, dy, x, y)
    return np.asarray((dphi_dx, dphi_dy))
# phi-field-spatial-gradients ends here
# [[file:../../org/pyretechnics.org::phi-east][phi-east]]
def calc_phi_east(phi, u_x, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y][x+1] - phi[y][x]
    if u_x >= 0:
        dphi_up = phi[y][x] - phi[y][x-1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y][x+2] - phi[y][x+1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x+1] - 0.5 * B * dphi_loc
# phi-east ends here
# [[file:../../org/pyretechnics.org::phi-west][phi-west]]
def calc_phi_west(phi, u_x, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y][x-1] - phi[y][x]
    if u_x >= 0:
        dphi_up = phi[y][x-2] - phi[y][x-1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x-1] - 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y][x] - phi[y][x+1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
# phi-west ends here
# [[file:../../org/pyretechnics.org::phi-north][phi-north]]
def calc_phi_north(phi, u_y, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y-1][x] - phi[y][x]
    if u_y >= 0:
        dphi_up = phi[y][x] - phi[y+1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y-2][x] - phi[y-1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y-1][x] - 0.5 * B * dphi_loc
# phi-north ends here
# [[file:../../org/pyretechnics.org::phi-south][phi-south]]
def calc_phi_south(phi, u_y, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y+1][x] - phi[y][x]
    if u_y >= 0:
        dphi_up = phi[y+2][x] - phi[y+1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y+1][x] - 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y][x] - phi[y-1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
# phi-south ends here
# [[file:../../org/pyretechnics.org::phi-time][phi-time]]
import numpy as np
from pyretechnics.burn_cells import burn_cell_toward_azimuth


# FIXME: stub
def identify_perimeter_cells(phi):
    """
    TODO: Add docstring
    """
    return [[0,0]]


def calc_phi_star(phi, u_x, u_y, dx, dy, dt, x, y):
    """
    Return an estimate for phi[y][x] at time (t + dt) given these inputs:
    - phi :: 2D float array of values in [-1,1]
    - u_x :: m/min
    - u_y :: m/min
    - dx  :: meters
    - dy  :: meters
    - dt  :: minutes
    - x   :: integer column index in phi
    - y   :: integer row index in phi
    """
    spread_vector   = np.asarray((u_x, u_y))
    gradient_vector = calc_phi_gradient(phi, u_x, u_y, dx, dy, x, y)
    dphi_dt         = np.dot(spread_vector, gradient_vector)
    return phi[y][x] - dphi_dt * dt


# TODO: Add exit conditions for a max_duration or no perimeter cells found and that should be the spread algorithm!
# TODO: Store the fire behavior values computed in the first perimeter_cells loop in output arrays.
def calc_phi_next_timestep(space_time_cubes, phi, dx, dy, dt, t):
    """
    TODO: Add docstring
    NOTE:
    - space_time_cubes and phi must have the same spatial resolution.
    - space_time_cubes must support temporal lookups in minutes.
    - dx is the cell width in meters.
    - dy is the cell height in meters.
    - dt is the timestep in minutes.
    - t is the start time in minutes.
    """
    perimeter_cells = identify_perimeter_cells(phi)

    phi_star = np.copy(phi)
    # TODO: Make this into a ufunc and apply directly to the array
    for [y, x] in perimeter_cells:
        # Calculate the spread vector normal to the fire front on the slope-tangential plane
        # FIXME: This only gives the wavelet spread rate, not the fire front spread rate. Use Val's functions instead.
        space_time_coordinate = (t, y, x)
        normal_vector         = calc_phi_normal_vector(phi, dx, dy, x, y)
        normal_azimuth        = calc_phi_normal_azimuth(normal_vector)
        fire_behavior         = burn_cell_toward_azimuth(space_time_cubes, space_time_coordinate, normal_azimuth)
        (u_x, u_y, u_z)       = fire_behavior["spread_rate"] * fire_behavior["spread_direction"]

        # Calculate the gradient of phi given this spread vector projected onto the horizontal plane
        # Update phi_star based on the dot product of the spread vector and the gradient vector
        phi_star[y][x] = calc_phi_star(phi, u_x, u_y, dx, dy, dt, x, y)

    perimeter_cells_star = identify_perimeter_cells(phi_star)

    phi_star_star = np.copy(phi_star)
    # TODO: Make this into a ufunc and apply directly to the array
    for [y, x] in perimeter_cells_star:
        # Calculate the spread vector normal to the fire front on the slope-tangential plane
        # FIXME: This only gives the wavelet spread rate, not the fire front spread rate. Use Val's functions instead.
        space_time_coordinate = (t + dt, y, x)
        normal_vector         = calc_phi_normal_vector(phi_star, dx, dy, x, y)
        normal_azimuth        = calc_phi_normal_azimuth(normal_vector)
        fire_behavior         = burn_cell_toward_azimuth(space_time_cubes, space_time_coordinate, normal_azimuth)
        (u_x, u_y, u_z)       = fire_behavior["spread_rate"] * fire_behavior["spread_direction"]

        # Calculate the gradient of phi_star given this spread vector projected onto the horizontal plane
        # Update phi_star_star based on the dot product of the spread vector and the gradient vector
        phi_star_star[y][x] = calc_phi_star(phi_star, u_x, u_y, dx, dy, dt, x, y)

    return (phi + phi_star_star) / 2.0
# phi-time ends here
# [[file:../../org/pyretechnics.org::calc-fireline-normal-behavior][calc-fireline-normal-behavior]]
from math import sqrt
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.surface_fire as sf
import pyretechnics.vector_utils as vu


# TODO: Move this to pyretechnics.vector_utils and use throughout the literate program
def calc_elevation_gradient(slope, aspect):
    """
    Returns the elevation gradient (dz_dx: rise/run, dz_dy: rise/run) given:
    - slope  :: rise/run
    - aspect :: degrees clockwise from North
    """
    return conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))


def calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient):
    """
    Return the gradient of phi projected onto the slope-tangential plane as a 3D (x,y,z) vector (in phi/m) given:
    - phi_gradient_xy    :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - elevation_gradient :: (dz_dx: m/m, dz_dy: m/m)
    """
    (dphi_dx, dphi_dy) = phi_gradient_xy
    phi_gradient_xyz   = np.asarray((dphi_dx, dphi_dy, 0.0))
    if vu.vector_magnitude(elevation_gradient) == 0.0:
        return phi_gradient_xyz
    else:
        slope_normal_vector = vu.get_slope_normal_vector(elevation_gradient) # (x,y,z) unit vector
        return phi_gradient_xyz - np.dot(phi_gradient_xyz, slope_normal_vector) * slope_normal_vector


# FIXME: Do I switch to cruz_passive_crown_fire_spread_rate() if the normal_spread_rate < critical_spread_rate?
#        Did I do this correctly in calc_crown_fire_behavior_in_direction?
def calc_fireline_normal_behavior(fire_behavior_max, phi_gradient):
    """
    Given these inputs:
    - fire_behavior_max  :: dictionary of max surface or crown fire behavior values
      - max_fire_type          :: "passive_crown" or "active_crown" (Required for crown fires only)
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_fireline_intensity :: kW/m
      - max_flame_length       :: m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
      - critical_spread_rate   :: m/min (Required for crown fires only)
    - phi_gradient       :: (dphi_dx: phi/m, dphi_dy: phi/m, dphi_dz: phi/m) 3D vector on the slope-tangential plane

    return a dictionary containing these keys:
    - dphi_dt            :: phi/min (on the slope-tangential plane)
    - fire_type          :: "unburned", "surface", "passive_crown", or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m

    Note: This function should work for surface or crown fires interchangeably.
    """
    #================================================================================================
    # Calculate the magnitude of the phi gradient
    #================================================================================================

    phi_magnitude = vu.vector_magnitude(phi_gradient) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burning
    #================================================================================================

    if (phi_magnitude == 0.0 or fire_behavior_max["max_spread_rate"] == 0.0):
        # This location is not on the fire perimeter and/or is not burning

        #================================================================================================
        # Set the spread direction to the phi gradient direction, max spread direction, upslope, or North
        #================================================================================================

        spread_direction = (np.asarray(phi_gradient) / phi_magnitude if phi_magnitude > 0.0
                            else fire_behavior_max["max_spread_direction"])

        #============================================================================================
        # Return zero surface/crown fire behavior
        #============================================================================================

        return {
            "dphi_dt"           : 0.0,
            "fire_type"         : "unburned",
            "spread_rate"       : 0.0,
            "spread_direction"  : spread_direction,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }

    else:
        # This location is on the fire perimeter and is burning

        #============================================================================================
        # Unpack the fire_behavior_max dictionary
        #============================================================================================

        heading_fire_type          = fire_behavior_max.get("max_fire_type", "surface")
        heading_spread_rate        = fire_behavior_max["max_spread_rate"]               # m/min
        heading_spread_direction   = fire_behavior_max["max_spread_direction"]          # (x,y,z) unit vector
        heading_spread_vector      = heading_spread_rate * heading_spread_direction     # (x,y,z) m/min vector
        heading_fireline_intensity = fire_behavior_max["max_fireline_intensity"]        # kW/m
        length_to_width_ratio      = fire_behavior_max["length_to_width_ratio"]         # unitless
        eccentricity               = fire_behavior_max["eccentricity"]                  # unitless
        critical_spread_rate       = fire_behavior_max.get("critical_spread_rate", 0.0) # m/min

        #============================================================================================
        # Calculate the backing and flanking fire spread rates
        #============================================================================================

        backing_adjustment   = (1.0 - eccentricity) / (1.0 + eccentricity)                                 # unitless
        backing_spread_rate  = heading_spread_rate * backing_adjustment                                    # m/min
        flanking_spread_rate = (heading_spread_rate + backing_spread_rate) / (2.0 * length_to_width_ratio) # m/min

        #============================================================================================
        # Calculate dphi/dt
        #============================================================================================

        A       = (heading_spread_rate - backing_spread_rate) / (2 * heading_spread_rate) # unitless
        B       = np.dot(heading_spread_vector, phi_gradient)                             # phi/min
        C       = flanking_spread_rate / heading_spread_rate                              # unitless
        D       = (heading_spread_rate * phi_magnitude) ** 2.0                            # (phi/min)^2
        E       = (length_to_width_ratio ** 2.0 - 1.0) * (B ** 2.0)                       # (phi/min)^2
        dphi_dt = -(A * B + C * sqrt(D + E))                                              # phi/min

        #============================================================================================
        # Calculate fire behavior normal to the fire perimeter
        #============================================================================================

        normal_spread_rate        = -dphi_dt / phi_magnitude                        # m/min
        normal_direction          = np.asarray(phi_gradient) / phi_magnitude        # (x,y,z) unit vector
        normal_adjustment         = normal_spread_rate / heading_spread_rate        # unitless
        normal_fireline_intensity = heading_fireline_intensity * normal_adjustment  # kW/m
        normal_flame_length       = sf.calc_flame_length(normal_fireline_intensity) # m
        normal_fire_type          = ("surface" if heading_fire_type == "surface"
                                     else "active_crown" if normal_spread_rate > critical_spread_rate
                                     else "passive_crown")

        #========================================================================================
        # Return the surface/crown fire behavior normal to the fire perimeter
        #========================================================================================

        return {
            "dphi_dt"           : dphi_dt,                   # phi/min
            "fire_type"         : normal_fire_type,          # surface, passive_crown, or active_crown
            "spread_rate"       : normal_spread_rate,        # m/min
            "spread_direction"  : normal_direction,          # (x,y,z) unit vector
            "fireline_intensity": normal_fireline_intensity, # kW/m
            "flame_length"      : normal_flame_length,       # m
        }
# calc-fireline-normal-behavior ends here
# [[file:../../org/pyretechnics.org::burn-cell-toward-phi-gradient][burn-cell-toward-phi-gradient]]
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.crown_fire as cf
import pyretechnics.fuel_models as fm
import pyretechnics.surface_fire as sf
import pyretechnics.vector_utils as vu


# TODO: Create a version of this function that runs efficiently over a space_time_region
def burn_cell_toward_phi_gradient(space_time_cubes, space_time_coordinate, phi_gradient_xy,
                                  use_wind_limit=True, max_length_to_width_ratio=None):
    """
    Given these inputs:
    - space_time_cubes          :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_10hr       :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_100hr      :: kg moisture/kg ovendry weight
      - fuel_moisture_live_herbaceous :: kg moisture/kg ovendry weight
      - fuel_moisture_live_woody      :: kg moisture/kg ovendry weight
      - foliar_moisture               :: kg moisture/kg ovendry weight
      - fuel_spread_adjustment        :: float >= 0.0 (Optional: defaults to 1.0)
      - weather_spread_adjustment     :: float >= 0.0 (Optional: defaults to 1.0)
    - space_time_coordinate     :: (t,y,x)
    - phi_gradient_xy           :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - use_wind_limit            :: boolean (Optional)
    - max_length_to_width_ratio :: float > 0.0 (Optional)

    return a dictionary with these fire behavior values for the space-time coordinate (t,y,x):
    - dphi_dt            :: phi/min (on the slope-tangential plane)
    - fire_type          :: "unburned", "surface", "passive_crown", or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector on the slope-tangential plane
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    #================================================================================================
    # Destructure the space_time_coordinate
    #================================================================================================

    (t, y, x) = space_time_coordinate

    #================================================================================================
    # Unpack the space_time_cubes dictionary
    #================================================================================================

    # Topography, Fuel Model, and Vegetation
    slope               = space_time_cubes["slope"].get(t,y,x)               # rise/run
    aspect              = space_time_cubes["aspect"].get(t,y,x)              # degrees clockwise from North
    fuel_model_number   = space_time_cubes["fuel_model"].get(t,y,x)          # integer index in fm.fuel_model_table
    canopy_cover        = space_time_cubes["canopy_cover"].get(t,y,x)        # 0-1
    canopy_height       = space_time_cubes["canopy_height"].get(t,y,x)       # m
    canopy_base_height  = space_time_cubes["canopy_base_height"].get(t,y,x)  # m
    canopy_bulk_density = space_time_cubes["canopy_bulk_density"].get(t,y,x) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                = space_time_cubes["wind_speed_10m"].get(t,y,x)                # km/hr
    upwind_direction              = space_time_cubes["upwind_direction"].get(t,y,x)              # degrees clockwise from North
    fuel_moisture_dead_1hr        = space_time_cubes["fuel_moisture_dead_1hr"].get(t,y,x)        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       = space_time_cubes["fuel_moisture_dead_10hr"].get(t,y,x)       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      = space_time_cubes["fuel_moisture_dead_100hr"].get(t,y,x)      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous = space_time_cubes["fuel_moisture_live_herbaceous"].get(t,y,x) # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      = space_time_cubes["fuel_moisture_live_woody"].get(t,y,x)      # kg moisture/kg ovendry weight
    foliar_moisture               = space_time_cubes["foliar_moisture"].get(t,y,x)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    = (space_time_cubes["fuel_spread_adjustment"].get(t,y,x)
                                 if "fuel_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    weather_spread_adjustment = (space_time_cubes["weather_spread_adjustment"].get(t,y,x)
                                 if "weather_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    spread_rate_adjustment    = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

    #================================================================================================
    # Calculate the elevation gradient
    #================================================================================================

    elevation_gradient = calc_elevation_gradient(slope, aspect)

    #============================================================================================
    # Project the horizontal phi gradient onto the slope-tangential plane
    #============================================================================================

    phi_gradient = calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient)

    #================================================================================================
    # Calculate the magnitude of the phi gradient
    #================================================================================================

    phi_magnitude = vu.vector_magnitude(phi_gradient) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burnable
    #================================================================================================

    fuel_model = fm.fuel_model_table.get(fuel_model_number)

    if not (phi_magnitude > 0.0 and fuel_model and fuel_model["burnable"]):
        # Cell is not on the fire perimeter and/or contains an unknown or non-burnable fuel model

        #================================================================================================
        # Set the spread direction to the phi gradient direction, upslope, or North
        #================================================================================================

        if phi_magnitude > 0.0:
            spread_direction = np.asarray(phi_gradient) / phi_magnitude
        elif slope > 0.0:
            slope_vector_3d  = vu.to_slope_plane(elevation_gradient, elevation_gradient)
            spread_direction = vu.as_unit_vector(slope_vector_3d)
        else:
            spread_direction = np.asarray((0,1,0)) # default: North

        #============================================================================================
        # Return zero surface fire behavior
        #============================================================================================

        return {
            "dphi_dt"           : 0.0,
            "fire_type"         : "unburned",
            "spread_rate"       : 0.0,
            "spread_direction"  : spread_direction,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }

    else:
        # Cell is on the fire perimeter and contains a burnable fuel model

        #============================================================================================
        # Compute derived parameters
        #============================================================================================

        fuel_moisture                = [fuel_moisture_dead_1hr,
                                        fuel_moisture_dead_10hr,
                                        fuel_moisture_dead_100hr,
                                        0.0, # fuel_moisture_dead_herbaceous
                                        fuel_moisture_live_herbaceous,
                                        fuel_moisture_live_woody]               # kg moisture/kg ovendry weight
        fuel_bed_depth               = fuel_model["delta"]                      # ft
        heat_of_combustion           = conv.Btu_lb_to_kJ_kg(fuel_model["h"][0]) # kJ/kg
        estimated_fine_fuel_moisture = fuel_moisture_dead_1hr                   # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        midflame_wind_speed = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,       # m/min
                                                          fuel_bed_depth,              # ft
                                                          conv.m_to_ft(canopy_height), # ft
                                                          canopy_cover)                # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        moisturized_fuel_model = fm.moisturize(fuel_model, fuel_moisture)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                          spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                             midflame_wind_speed,
                                                             upwind_direction,
                                                             slope,
                                                             aspect,
                                                             use_wind_limit)

        #============================================================================================
        # Calculate surface fire behavior normal to the fire perimeter
        #============================================================================================

        surface_fire_normal = calc_fireline_normal_behavior(surface_fire_max, phi_gradient)

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(surface_fire_normal["fireline_intensity"],
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            crown_fire_max = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                             canopy_bulk_density, heat_of_combustion,
                                                             estimated_fine_fuel_moisture,
                                                             wind_speed_10m, upwind_direction,
                                                             slope, aspect, max_length_to_width_ratio)

            #========================================================================================
            # Calculate crown fire behavior normal to the fire perimeter
            #========================================================================================

            crown_fire_normal = calc_fireline_normal_behavior(crown_fire_max, phi_gradient)

            #========================================================================================
            # Calculate combined fire behavior normal to the fire perimeter
            #========================================================================================

            combined_fire_normal = cf.calc_combined_fire_behavior(surface_fire_normal, crown_fire_normal)
            surface_dphi_dt      = surface_fire_normal["dphi_dt"]
            crown_dphi_dt        = crown_fire_normal["dphi_dt"]
            combined_dphi_dt     = surface_dphi_dt if abs(surface_dphi_dt) > abs(crown_dphi_dt) else crown_dphi_dt
            combined_fire_normal["dphi_dt"] = combined_dphi_dt

            #========================================================================================
            # Return the combined fire behavior normal to the fire perimeter
            #========================================================================================

            return combined_fire_normal

        else:

            #========================================================================================
            # Return the surface fire behavior normal to the fire perimeter
            #========================================================================================

            return surface_fire_normal
# burn-cell-toward-phi-gradient ends here
# [[file:../../org/pyretechnics.org::identify-perimeter-cells-from-phi-field][identify-perimeter-cells-from-phi-field]]
def opposite_phi_signs(phi_matrix, y1, x1, y2, x2):
    """
    TODO: Add docstring
    """
    return phi_matrix[y1, x1] * phi_matrix[y2, x2] < 0.0


def identify_perimeter_cells(phi_matrix):
    """
    TODO: Add docstring
    """
    frontier_cells = []

    (rows, cols) = phi_matrix.shape

    # Scan interior cells
    for y in range(1, rows-1):
        for x in range(1, cols-1):
            # Compare (north, south) and (east, west) neighboring cell pairs for opposite phi signs
            if opposite_phi_signs(phi_matrix, y-1, x, y+1, x) or opposite_phi_signs(phi_matrix, y, x+1, y, x-1):
                frontier_cells.append((y, x))

    # Scan northern edge (non-corner) cells
    y = 0
    for x in range(1, cols-1):
        # Compare (here, south) and (east, west) neighboring cell pairs for opposite phi signs
        if opposite_phi_signs(phi_matrix, y, x, y+1, x) or opposite_phi_signs(phi_matrix, y, x+1, y, x-1):
            frontier_cells.append((y, x))

    # Scan southern edge (non-corner) cells
    y = rows-1
    for x in range(1, cols-1):
        # Compare (here, north) and (east, west) neighboring cell pairs for opposite phi signs
        if opposite_phi_signs(phi_matrix, y, x, y-1, x) or opposite_phi_signs(phi_matrix, y, x+1, y, x-1):
            frontier_cells.append((y, x))

    # Scan eastern edge (non-corner) cells
    x = cols-1
    for y in range(1, rows-1):
        # Compare (north, south) and (here, west) neighboring cell pairs for opposite phi signs
        if opposite_phi_signs(phi_matrix, y-1, x, y+1, x) or opposite_phi_signs(phi_matrix, y, x, y, x-1):
            frontier_cells.append((y, x))

    # Scan western edge (non-corner) cells
    x = 0
    for y in range(1, rows-1):
        # Compare (north, south) and (here, east) neighboring cell pairs for opposite phi signs
        if opposite_phi_signs(phi_matrix, y-1, x, y+1, x) or opposite_phi_signs(phi_matrix, y, x, y, x+1):
            frontier_cells.append((y, x))

    # Scan northwestern corner cell
    y = 0
    x = 0
    # Compare (here, south) and (here, east) neighboring cell pairs for opposite phi signs
    if opposite_phi_signs(phi_matrix, y, x, y+1, x) or opposite_phi_signs(phi_matrix, y, x, y, x+1):
        frontier_cells.append((y, x))

    # Scan northeastern corner cell
    y = 0
    x = cols-1
    # Compare (here, south) and (here, west) neighboring cell pairs for opposite phi signs
    if opposite_phi_signs(phi_matrix, y, x, y+1, x) or opposite_phi_signs(phi_matrix, y, x, y, x-1):
        frontier_cells.append((y, x))

    # Scan southwestern corner cell
    y = rows-1
    x = 0
    # Compare (here, north) and (here, east) neighboring cell pairs for opposite phi signs
    if opposite_phi_signs(phi_matrix, y, x, y-1, x) or opposite_phi_signs(phi_matrix, y, x, y, x+1):
        frontier_cells.append((y, x))

    # Scan southeastern corner cell
    y = rows-1
    x = cols-1
    # Compare (here, north) and (here, west) neighboring cell pairs for opposite phi signs
    if opposite_phi_signs(phi_matrix, y, x, y-1, x) or opposite_phi_signs(phi_matrix, y, x, y, x-1):
        frontier_cells.append((y, x))

    # Return the list of frontier cells
    return frontier_cells
# identify-perimeter-cells-from-phi-field ends here
# [[file:../../org/pyretechnics.org::spread-phi-field][spread-phi-field]]
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.vector_utils as vu


# TODO: Move to pyretechnics.conversion
fire_type_codes = {
    "unburned"      : 0,
    "surface"       : 1,
    "passive_crown" : 2,
    "active_crown"  : 3,
}


# TODO: Move to pyretechnics.vector_utils
def spread_direction_vector_to_angle(vector_3d):
    """
    TODO: Add docstring
    """
    if vector_3d:
        (x, y)       = vu.to_horizontal_plane(vector_3d)
        (r, azimuth) = conv.cartesian_to_azimuthal(x, y)
        return azimuth
    else:
        return 0.0 # default: North


def spread_fire_one_timestep(space_time_cubes, output_matrices, cell_width, cell_height, start_time,
                             use_wind_limit=True, max_length_to_width_ratio=None, max_cells_per_timestep=1):
    """
    TODO: Add docstring
    NOTE:
    - space_time_cubes and phi must have the same spatial resolution and extent.
    - space_time_cubes must support temporal lookups in minutes.
    - cell_width is in meters.
    - cell_height is in meters.
    - dt is the timestep in minutes.
    - start_time is the start time in minutes.
    """
    # Unpack output matrices
    phi_matrix                = output_matrices["phi"]
    fire_type_matrix          = output_matrices["fire_type"]
    spread_rate_matrix        = output_matrices["spread_rate"]
    spread_direction_matrix   = output_matrices["spread_direction"]
    fireline_intensity_matrix = output_matrices["fireline_intensity"]
    flame_length_matrix       = output_matrices["flame_length"]
    time_of_arrival_matrix    = output_matrices["time_of_arrival"]

    # TODO: Compute max_timestep based on the temporal resolution of the space_time_cubes
    max_timestep = 1.0

    # Initialize max spread rates in the x and y dimensions to 0.0
    max_spread_rate_x = 0.0
    max_spread_rate_y = 0.0

    # Create an empty dictionary to store intermediate fire behavior values per cell
    fire_behavior_dict = {}

    # Identify perimeter cells
    perimeter_cells = identify_perimeter_cells(phi_matrix)

    # Make a copy of phi_matrix as phi_star_matrix
    phi_star_matrix = np.copy(phi_matrix)

    # Compute fire behavior values at time (start_time) and identify the max spread rates in the x and y dimensions
    for cell_index in perimeter_cells:
        # Unpack cell_index
        (y, x) = cell_index

        # Calculate phi gradient on the horizontal plane
        phi_gradient_xy = calc_phi_gradient_approx(phi_matrix, cell_width, cell_height, x, y)

        # Calculate the fire behavior normal to the fire front on the slope-tangential plane
        fire_behavior = burn_cell_toward_phi_gradient(space_time_cubes, (start_time, y, x), phi_gradient_xy,
                                                      use_wind_limit, max_length_to_width_ratio)

        # Keep a running tally of the max spread rates in the x and y dimensions
        (spread_rate_x, spread_rate_y, _) = fire_behavior["spread_rate"] * fire_behavior["spread_direction"]
        max_spread_rate_x = max(max_spread_rate_x, abs(spread_rate_x))
        max_spread_rate_y = max(max_spread_rate_y, abs(spread_rate_y))

        # Integrate the Superbee flux limited phi gradient to make dphi_dt numerically stable
        phi_magnitude = vu.vector_magnitude(phi_gradient_xy)
        if phi_magnitude > 0.0:
            phi_gradient_xy_limited = calc_phi_gradient(phi_matrix, *phi_gradient_xy, cell_width, cell_height, x, y)
            fire_behavior["dphi_dt"] *= np.dot(phi_gradient_xy, phi_gradient_xy_limited) / phi_magnitude

        # Store fire behavior values for later use
        fire_behavior_dict[cell_index] = fire_behavior

    # Calculate timestep using the CFL condition
    # TODO: Incorporate max_cells_per_timestep
    if max_spread_rate_x == 0.0:
        if max_spread_rate_y == 0.0:
            dt = max_timestep
        else:
            dt = cell_height / max_spread_rate_y
    else:
        if max_spread_rate_y == 0.0:
            dt = cell_width / max_spread_rate_x
        else:
            dt = min(cell_width / max_spread_rate_x,
                     cell_height / max_spread_rate_y)

    # Update the perimeter cell values in phi_star_matrix
    for cell_index in perimeter_cells:
        (y, x) = cell_index
        phi_star_matrix[y][x] += fire_behavior_dict[cell_index]["dphi_dt"] * dt

    # Compute fire behavior values at time (start_time + dt) and update the output matrices
    for cell_index in perimeter_cells:
        # Unpack cell_index
        (y, x) = cell_index

        # Calculate phi gradient on the horizontal plane
        phi_gradient_xy_star = calc_phi_gradient_approx(phi_star_matrix, cell_width, cell_height, x, y)

        # Calculate the fire behavior normal to the fire front on the slope-tangential plane
        fire_behavior_star = burn_cell_toward_phi_gradient(space_time_cubes, (start_time + dt, y, x),
                                                           phi_gradient_xy_star, use_wind_limit,
                                                           max_length_to_width_ratio)

        # Integrate the Superbee flux limited phi gradient to make dphi_dt numerically stable
        phi_magnitude = vu.vector_magnitude(phi_gradient_xy_star)
        if phi_magnitude > 0.0:
            phi_gradient_xy_star_limited = calc_phi_gradient(phi_star_matrix, *phi_gradient_xy_star,
                                                             cell_width, cell_height, x, y)
            fire_behavior_star["dphi_dt"] *= np.dot(phi_gradient_xy_star, phi_gradient_xy_star_limited) / phi_magnitude

        # Calculate the new phi value at time (start_time + dt) as phi_next
        fire_behavior     = fire_behavior_dict[cell_index]
        dphi_dt_estimate1 = fire_behavior["dphi_dt"]
        dphi_dt_estimate2 = fire_behavior_star["dphi_dt"]
        dphi_dt_average   = (dphi_dt_estimate1 + dphi_dt_estimate2) / 2.0
        phi               = phi_matrix[y][x]
        phi_next          = phi + dphi_dt_average * dt

        # Update the perimeter cell values in phi_matrix
        phi_matrix[y][x] = phi_next

        # Record fire behavior values in their respective output matrices for cells that are burned in this timestep
        # NOTE: This records the fire behavior values at time (start_time) and not at the time of arrival.
        if phi > 0.0 and phi_next <= 0.0:
            fire_type_matrix[y][x]          = fire_type_codes[fire_behavior["fire_type"]]
            spread_rate_matrix[y][x]        = fire_behavior["spread_rate"]
            spread_direction_matrix[y][x]   = spread_direction_vector_to_angle(fire_behavior["spread_direction"])
            fireline_intensity_matrix[y][x] = fire_behavior["fireline_intensity"]
            flame_length_matrix[y][x]       = fire_behavior["flame_length"]
            time_of_arrival_matrix[y][x]    = start_time + dt * phi / (phi - phi_next)

    return {
        "simulation_time": start_time + dt,
        "output_matrices": output_matrices,
    }


# FIXME: stub
# TODO: Add exit conditions for a max_duration or no perimeter cells found and that should be the spread algorithm!
def spread_fire_with_phi_field(space_time_cubes, output_matrices, cell_width, cell_height, start_time,
                               use_wind_limit=True, max_length_to_width_ratio=None, max_cells_per_timestep=1):
    """
    Given these inputs:
    - space_time_cubes          :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - wind_speed_10m                :: km/hr
      - upwind_direction              :: degrees clockwise from North
      - fuel_moisture_dead_1hr        :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_10hr       :: kg moisture/kg ovendry weight
      - fuel_moisture_dead_100hr      :: kg moisture/kg ovendry weight
      - fuel_moisture_live_herbaceous :: kg moisture/kg ovendry weight
      - fuel_moisture_live_woody      :: kg moisture/kg ovendry weight
      - foliar_moisture               :: kg moisture/kg ovendry weight
      - fuel_spread_adjustment        :: float >= 0.0 (Optional: defaults to 1.0)
      - weather_spread_adjustment     :: float >= 0.0 (Optional: defaults to 1.0)
    - phi                       :: 2D float array of values in [-1,1]
    - cell_width                :: meters
    - cell_height               :: meters
    - start_time                :: minutes
    - use_wind_limit            :: boolean (Optional)
    - max_length_to_width_ratio :: float > 0.0 (Optional)

    return a dictionary with these fire behavior values for the space-time region:
    - phi_final          :: 2D float array of values in [-1,1]
    - fire_type          :: 2D byte array (0=unburned, 1=surface, 2=passive_crown, 3=active_crown)
    - spread_rate        :: 2D float array (m/min)
    - spread_direction   :: 2D float array (degrees clockwise from North)
    - fireline_intensity :: 2D float array (kW/m)
    - flame_length       :: 2D float array (m)
    - time_of_arrival    :: 2D float array (min)
    """
    # Sketch of Fire Spread Progression Algorithm:
    # ==============================================================================================
    # 2. Exit if start_time + dt > max_stop_time.
    # 3. Create phi_star and phi_next 2D arrays.
    # 4. Identify perimeter cells with a buffer from phi 2D array.
    # 5. Exit if no perimeter cells are found.
    # 6. For each perimeter cell [y,x]:
    #    3.1. phi_gradient_xy = calc_phi_gradient_approx(phi, cell_width, cell_height, x, y)
    #    3.2. if vector_magnitude(phi_gradient_xy) == 0.0:
    #         3.2.1. Skip this cell since it is not on the perimeter. Return to 3.1. for the next perimeter cell.
    #         else:
    #         3.2.2  Look up slope and aspect from the space_time_cubes dictionary.
    #         3.2.3. elevation_gradient = calc_elevation_gradient(slope, aspect)
    #         3.2.4. phi_gradient_on_slope = calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient)
    #         3.2.5. (dphi_dt, fire_behavior) = burn_cell_toward_phi_gradient(space_time_cubes, space_time_coordinate,
    #                                                                         phi_gradient, use_wind_limit,
    #                                                                         max_length_to_width_ratio)
    #    3.3. Calculate dt using the CFL condition.
    #    3.4. Set phi_star[y][x] = phi[y][x] + dphi_dt * dt
    # 6. Jump to step 4 and set the time to (start_time + dt). Also replace phi -> phi_star and phi_star -> phi_next.
    # 7. Set the phi 2D array to be the average of phi and phi_next.
    # 8. Jump to step 2 and repeat for the next timestep.

    results = spread_fire_one_timestep(space_time_cubes, output_matrices, cell_width, cell_height, start_time,
                                       use_wind_limit=True, max_length_to_width_ratio=None, max_cells_per_timestep=1)

    new_time                = results["simulation_time"]
    updated_output_matrices = results["output_matrices"]

    return None
# spread-phi-field ends here
# [[file:../../org/pyretechnics.org::*Spread Phi Field][Spread Phi Field:3]]
# C_max will usually be 1
dt = min(C_max * cell_size / max(max(x_spread_rates), max(y_spread_rates)), simulation_dtmax)
# Spread Phi Field:3 ends here
