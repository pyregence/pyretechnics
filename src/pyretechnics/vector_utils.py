# [[file:../../org/pyretechnics.org::vector-utilities][vector-utilities]]
from math import radians, sin, cos, sqrt
import numpy as np
import pyretechnics.conversion as conv


def vector_magnitude(vector):
    return sqrt(np.dot(vector, vector))


def as_unit_vector(vector):
    return np.asarray(vector) / vector_magnitude(vector)


def to_slope_plane(vector_2d, elevation_gradient):
    return np.asarray((
        vector_2d[0],
        vector_2d[1],
        np.dot(vector_2d, elevation_gradient)
    ))


def to_horizontal_plane(vector_3d):
    return np.asarray(vector_3d[0:2])


def get_slope_normal_vector(elevation_gradient):
    slope_normal_vector = (-elevation_gradient[0], -elevation_gradient[1], 1)
    return as_unit_vector(slope_normal_vector)


def rotate_on_sloped_plane(vector, theta, slope, aspect):
    """Rotate a 3D vector <x,y,z> theta degrees clockwise on the plane defined by the slope and aspect."""
    # Calculate the slope normal vector from the slope and aspect
    elevation_gradient  = conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))
    slope_normal_vector = get_slope_normal_vector(elevation_gradient)
    # Rotate theta degrees clockwise around the slope_normal_vector
    theta_rad = radians(theta)
    return cos(theta_rad) * vector + np.matmul(sin(theta_rad) * slope_normal_vector, vector)
# vector-utilities ends here
