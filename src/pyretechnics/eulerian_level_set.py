# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients][phi-field-spatial-gradients]]
def calc_dphi_dx_approx(phi, x, y, dx):
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given the cell width dx.

    NOTE: The origin cell (x=0,y=0) is located in the upper left corner
          of the grid in Python arrays. Thus, as x increases, we move
          to the east, and as y increases, we move to the south.
    """
    return (phi[y][x+1] - phi[y][x-1]) / (2.0 * dx)


def calc_dphi_dy_approx(phi, x, y, dy):
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given the cell height dy.

    NOTE: The origin cell (x=0,y=0) is located in the upper left corner
          of the grid in Python arrays. Thus, as x increases, we move
          to the east, and as y increases, we move to the south.
    """
    return (phi[y-1][x] - phi[y+1][x]) / (2.0 * dy)
# phi-field-spatial-gradients ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector][phi-field-normal-vector]]
from math import sqrt

def calc_phi_gradient_magnitude(dphi_dx, dphi_dy):
    """Calculate the magnitude of the gradient of the phi field."""
    return sqrt(dphi_dx ** 2 + dphi_dy ** 2)


def calc_phi_normal_vector(phi, dx, dy, x, y):
    """
    Calculate the phi field normal vector in the x and y dimensions.

    n_x: eastward component of the unit normal vector
    n_y: northward component of the unit normal vector
    """
    dphi_dx = calc_dphi_dx_approx(phi, x, y, dx)
    dphi_dy = calc_dphi_dy_approx(phi, x, y, dy)
    if dphi_dx == 0.0 and dphi_dy == 0.0:
        return {
            "n_x": 0.0,
            "n_y": 0.0,
        }
    else:
        phi_gradient_magnitude = calc_phi_gradient_magnitude(dphi_dx, dphi_dy)
        return {
            "n_x": dphi_dx / phi_gradient_magnitude,
            "n_y": dphi_dy / phi_gradient_magnitude,
        }
# phi-field-normal-vector ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector-angle][phi-field-normal-vector-angle]]
from math import atan2, pi

def calc_phi_normal_vector_angle(phi_normal_vector):
    """
    Calculate the angle (measured in radians clockwise from North)
    to which the phi field's normal vector points.

    NOTE: The origin cell (x=0,y=0) is located in the upper left corner
          of the grid in Python arrays. Thus, as x increases, we move
          to the east, and as y increases, we move to the south.
    """
    n_x = phi_normal_vector["n_x"]
    n_y = phi_normal_vector["n_y"]
    if n_x >= 0 and n_y >= 0:
        return 1/2 * pi - atan2(n_y, n_x)
    elif n_x < 0 and n_y >= 0:
        return 3/2 * pi + atan2(n_y, abs(n_x))
    elif n_x < 0 and n_y < 0:
        return 3/2 * pi - atan2(n_y, n_x)
    elif n_x >= 0 and n_y < 0:
        return 1/2 * pi + atan2(abs(n_y), n_x)
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
# [[file:../../org/pyretechnics.org::dphi_partial_derivatives][dphi_partial_derivatives]]
def calc_dphi_dx(phi_east, phi_west, dx):
    """
    TODO: Add docstring
    """
    return (phi_east - phi_west) / dx

def calc_dphi_dy(phi_north, phi_south, dy):
    """
    TODO: Add docstring
    """
    return (phi_north - phi_south) / dy
# dphi_partial_derivatives ends here
# [[file:../../org/pyretechnics.org::phi_east][phi_east]]
def calc_phi_east(phi, u_x, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y][x+1] - phi[y][x]
    if u_x[y][x] >= 0:
        dphi_up = phi[y][x] - phi[y][x-1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y][x+2] - phi[y][x+1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x+1] - 0.5 * B * dphi_loc
# phi_east ends here
# [[file:../../org/pyretechnics.org::phi_west][phi_west]]
def calc_phi_west(phi, u_x, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y][x-1] - phi[y][x]
    if u_x[y][x] >= 0:
        dphi_up = phi[y][x-2] - phi[y][x-1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x-1] - 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y][x] - phi[y][x+1]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
# phi_west ends here
# [[file:../../org/pyretechnics.org::phi_north][phi_north]]
def calc_phi_north(phi, u_y, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y-1][x] - phi[y][x]
    if u_y[y][x] >= 0:
        dphi_up = phi[y][x] - phi[y+1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y-2][x] - phi[y-1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y-1][x] - 0.5 * B * dphi_loc
# phi_north ends here
# [[file:../../org/pyretechnics.org::phi_south][phi_south]]
def calc_phi_south(phi, u_y, x, y):
    """
    TODO: Add docstring
    """
    dphi_loc = phi[y+1][x] - phi[y][x]
    if u_y[y][x] >= 0:
        dphi_up = phi[y+2][x] - phi[y+1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y+1][x] - 0.5 * B * dphi_loc
    else:
        dphi_up = phi[y][x] - phi[y-1][x]
        B = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
# phi_south ends here
# [[file:../../org/pyretechnics.org::phi_time][phi_time]]
def calc_phi_time(dx, dy, dt):
    """
    TODO: Add docstring (computes phi_{t+dt})
    """
    # FIXME: Figure out how to compute all the new terms here:
    phi_t = 0
    U_x = 0
    phi_east_t = 0
    phi_west_t = 0
    U_y = 0
    phi_north_t = 0
    phi_south_t = 0
    phi_star = 0
    U_x_star = 0
    phi_east_star = 0
    phi_west_star = 0
    U_y_star = 0
    phi_north_star = 0
    phi_south_star = 0
    phi_star = phi_t - dt * (U_x * (phi_east_t - phi_west_t) / dx + U_y * (phi_north_t - phi_south_t) / dy)
    return 0.5 * phi_t + 0.5 * (phi_star - dt * (U_x_star * (phi_east_star - phi_west_star) / dx + U_y_star * (phi_north_star - phi_south_star) / dy))
# phi_time ends here
