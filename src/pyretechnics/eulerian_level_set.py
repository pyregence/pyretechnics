# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients][phi-field-spatial-gradients]]
def calc_dphi_dx(phi, x, y, dx):
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given the cell width dx.

    NOTE: The origin cell (x=0,y=0) is located in the upper left corner
          of the grid in Python arrays. Thus, as x increases, we move
          to the east, and as y increases, we move to the south.
    """
    return (phi[y][x+1] - phi[y][x-1]) / (2.0 * dx)


def calc_dphi_dy(phi, x, y, dy):
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
    dphi_dx = calc_dphi_dx(phi, x, y, dx)
    dphi_dy = calc_dphi_dy(phi, x, y, dy)
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
