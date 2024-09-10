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
