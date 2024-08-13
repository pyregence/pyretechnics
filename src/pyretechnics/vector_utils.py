# [[file:../../org/pyretechnics.org::vector-utilities][vector-utilities]]
from math import degrees, radians, sin, cos, acos, atan2, sqrt
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


def rotate_x(vector, theta):
    """Rotate a 3D vector <x,y,z> theta degrees clockwise around the x axis."""
    theta_rad       = radians(theta)
    rotation_matrix = [[1.0,            0.0,             0.0],
                       [0.0, cos(theta_rad), -sin(theta_rad)],
                       [0.0, sin(theta_rad),  cos(theta_rad)]]
    return np.dot(rotation_matrix, vector)


def rotate_y(vector, theta):
    """Rotate a 3D vector <x,y,z> theta degrees clockwise around the y axis."""
    theta_rad       = radians(theta)
    rotation_matrix = [[ cos(theta_rad), 0.0, sin(theta_rad)],
                       [            0.0, 1.0,            0.0],
                       [-sin(theta_rad), 0.0, cos(theta_rad)]]
    return np.dot(rotation_matrix, vector)


def rotate_z(vector, theta):
    """Rotate a 3D vector <x,y,z> theta degrees clockwise around the z axis."""
    theta_rad       = radians(theta)
    rotation_matrix = [[cos(theta_rad), -sin(theta_rad), 0.0],
                       [sin(theta_rad),  cos(theta_rad), 0.0],
                       [           0.0,             0.0, 1.0]]
    return np.dot(rotation_matrix, vector)


def rotate_on_sloped_plane(vector, theta, slope, aspect):
    """Rotate a 3D vector <x,y,z> theta degrees clockwise on the plane defined by the slope and aspect."""
    # Calculate the slope normal vector from the slope and aspect
    elevation_gradient  = conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))
    slope_normal_vector = get_slope_normal_vector(elevation_gradient)
    (n_Sx, n_Sy, n_Sz)  = slope_normal_vector
    # Rotate around z axis to align n_S with the y-z plane
    z_angle = degrees(atan2(n_Sx, n_Sy))
    vector1 = rotate_z(vector, z_angle)
    # Rotate around x axis to align n_S with the z axis
    x_angle = degrees(acos(n_Sz / vector_magnitude(slope_normal_vector)))
    vector2 = rotate_x(vector1, x_angle)
    # Now the sloped plane is the x-y plane!
    # Rotate around z axis (a.k.a. n_S) to turn the input vector theta degrees clockwise
    vector3 = rotate_z(vector2, theta)
    # Reverse the earlier x axis rotation
    vector4 = rotate_x(vector3, -x_angle)
    # Reverse the initial z axis rotation
    vector5 = rotate_z(vector4, -z_angle)
    return vector5
# vector-utilities ends here
