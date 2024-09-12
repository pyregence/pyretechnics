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
    if u_x > 0:
        dphi_up = phi[y][x] - phi[y][x-1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    elif u_x < 0:
        dphi_up = phi[y][x+2] - phi[y][x+1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x+1] - 0.5 * B * dphi_loc
    else:
        # FIXME: What is the correct value for this case? Update the equation above accordingly.
        return 0.0
# phi-east ends here
# [[file:../../org/pyretechnics.org::phi-west][phi-west]]
def calc_phi_west(phi, u_x, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y][x-1] - phi[y][x]
    if u_x > 0:
        dphi_up = phi[y][x-2] - phi[y][x-1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x-1] - 0.5 * B * dphi_loc
    elif u_x < 0:
        dphi_up = phi[y][x] - phi[y][x+1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        # FIXME: What is the correct value for this case? Update the equation above accordingly.
        return 0.0
# phi-west ends here
# [[file:../../org/pyretechnics.org::phi-north][phi-north]]
def calc_phi_north(phi, u_y, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y-1][x] - phi[y][x]
    if u_y > 0:
        dphi_up = phi[y][x] - phi[y+1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    elif u_y < 0:
        dphi_up = phi[y-2][x] - phi[y-1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y-1][x] - 0.5 * B * dphi_loc
    else:
        # FIXME: What is the correct value for this case? Update the equation above accordingly.
        return 0.0
# phi-north ends here
# [[file:../../org/pyretechnics.org::phi-south][phi-south]]
def calc_phi_south(phi, u_y, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y+1][x] - phi[y][x]
    if u_y > 0:
        dphi_up = phi[y+2][x] - phi[y+1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y+1][x] - 0.5 * B * dphi_loc
    elif u_y < 0:
        dphi_up = phi[y][x] - phi[y-1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        # FIXME: What is the correct value for this case? Update the equation above accordingly.
        return 0.0
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
# [[file:../../org/pyretechnics.org::vw-elliptical-solver][vw-elliptical-solver]]
from math import sqrt
import numpy as np
from pyretechnics.surface_fire import calc_flame_length
from pyretechnics.vector_utils import vector_magnitude


# FIXME: Do I switch to cruz_passive_crown_fire_spread_rate() if the normal_spread_rate < critical_spread_rate?
#        Did I do this correctly in calc_crown_fire_behavior_in_direction?
def calc_fireline_normal_behavior(fire_behavior_max, elevation_gradient, phi_gradient):
    """
    Given these inputs:
    - fire_behavior_max  :: dictionary of max surface or crown fire behavior values
      - max_fire_type          :: "passive_crown" or "active_crown" (Required for crown fires only)
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_spread_vector      :: (x: m/min, y: m/min, z: m/min)
      - max_fireline_intensity :: kW/m
      - max_flame_length       :: m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
      - critical_spread_rate   :: m/min (Required for crown fires only)
    - elevation_gradient :: (dz_dx: m/m, dz_dy: m/m)
    - phi_gradient       :: (dphi_dx: phi/m, dphi_dy: phi/m, dphi_dz: phi/m) 3D vector on the slope-tangential plane

    return a dictionary containing these keys:
    - dphi_dt            :: phi/min
    - fire_type          :: "unburned", "surface", "passive_crown", or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m

    Note: This function should work for surface or crown fires interchangeably.
    """
    phi_gradient  = np.asarray(phi_gradient)       # (x: phi/m, y: phi/m, z: phi/m)
    phi_magnitude = vector_magnitude(phi_gradient) # phi/m
    if phi_magnitude == 0.0:
        # This location is not on the fire perimeter.
        return {
            "dphi_dt"           : 0.0,
            "fire_type"         : "unburned",
            "spread_rate"       : 0.0,
            "spread_direction"  : None,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }
    else:
        # This location is on the fire perimeter.
        heading_spread_rate         = fire_behavior_max["max_spread_rate"]         # m/min
        heading_spread_vector       = fire_behavior_max["max_spread_vector"]       # (x: m/min, y: m/min, z: m/min)
        heading_spread_rate_x       = heading_spread_vector[0]                     # m/min
        heading_spread_rate_y       = heading_spread_vector[1]                     # m/min
        length_to_width_ratio       = fire_behavior_max["length_to_width_ratio"]   # unitless
        eccentricity                = fire_behavior_max["eccentricity"]            # unitless
        backing_adjustment          = (1.0 - eccentricity) / (1.0 + eccentricity)  # unitless
        backing_spread_rate         = heading_spread_rate * backing_adjustment     # m/min
        flanking_spread_rate        = ((heading_spread_rate + backing_spread_rate)
                                       / (2.0 * length_to_width_ratio))            # m/min
        heading_fireline_intensity  = fire_behavior_max["max_fireline_intensity"]  # kW/m
        heading_fire_type           = fire_behavior_max.get("max_fire_type", "surface")
        critical_spread_rate        = fire_behavior_max.get("critical_spread_rate", 0.0) # m/min
        (dz_dx, dz_dy)              = elevation_gradient                           # (x: m/m, y: m/m)
        (dphi_dx, dphi_dy, dphi_dz) = phi_gradient                                 # (x: phi/m, y: phi/m, z: phi/m)
        A = (heading_spread_rate - backing_spread_rate) / (2 * heading_spread_rate)              # unitless
        B = heading_spread_rate_x * dphi_dx + heading_spread_rate_y * dphi_dy                    # phi/min
        C = (flanking_spread_rate / heading_spread_rate)                                         # unitless
        D = (dphi_dx * dz_dx + dphi_dy * dz_dy) ** 2.0 / (1 + dz_dx ** 2.0 + dz_dy ** 2.0)       # (phi/m)^2
        E = dphi_dx ** 2.0 + dphi_dy ** 2.0 - D                                                  # (phi/m)^2
        F = (heading_spread_rate ** 2.0) * E + (length_to_width_ratio ** 2.0 - 1.0) * (B ** 2.0) # (phi/min)^2
        dphi_dt                   = -(A * B + C * sqrt(F))                                       # phi/min
        normal_spread_rate        = -dphi_dt / sqrt(E)                                           # m/min
        normal_direction          = phi_gradient / phi_magnitude                                 # (x,y,z) unit vector
        normal_adjustment         = normal_spread_rate / heading_spread_rate                     # unitless
        normal_fireline_intensity = heading_fireline_intensity * normal_adjustment               # kW/m
        normal_flame_length       = calc_flame_length(normal_fireline_intensity)                 # m
        normal_fire_type          = ("unburned" if heading_spread_rate == 0.0
                                     else "surface" if heading_fire_type == "surface"
                                     else "active_crown" if normal_spread_rate > critical_spread_rate
                                     else "passive_crown")
        return {
            "dphi_dt"           : dphi_dt,                   # phi/min
            "fire_type"         : normal_fire_type,          # surface, passive_crown, or active_crown
            "spread_rate"       : normal_spread_rate,        # m/min
            "spread_direction"  : normal_direction,          # (x,y,z) unit vector
            "fireline_intensity": normal_fireline_intensity, # kW/m
            "flame_length"      : normal_flame_length,       # m
        }
# vw-elliptical-solver ends here
