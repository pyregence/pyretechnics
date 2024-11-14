# [[file:../../org/pyretechnics.org::vector-utilities][vector-utilities]]
from math import radians
import cython
import cython as cy
import numpy as np
import pyretechnics.conversion as conv


# TODO: Fix error with importing pyretechnics.types
if cy.compiled:
    from cython.cimports.pyretechnics.math import sqrt, sin, cos
    from cython.cimports.pyretechnics.types import pyidx, vec_xy, vec_xyz
else:
    from math import sqrt, sin, cos
    from pyretechnics.types import pyidx, vec_xy, vec_xyz


@cy.profile(False)
@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def vector_magnitude(vector: cy.double[:]) -> cy.double:
    i  : pyidx
    x  : cy.double
    acc: cy.double = 0.0
    for i in range(len(vector)):
        x    = vector[i]
        acc += x * x
    return sqrt(acc)


# TODO: result uninitialized warning
@cy.profile(False)
@cy.ccall
def dot_2d(vector1: vec_xy, vector2: vec_xy) -> cy.float:
    return vector1[0] * vector2[0] + vector1[1] * vector2[1]


@cy.profile(False)
@cy.ccall
def dot_3d(vector1: vec_xyz, vector2: vec_xyz) -> cy.float:
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]


# TODO: result uninitialized warning
@cy.profile(False)
@cy.ccall
def vector_magnitude_2d(vector: vec_xy) -> cy.float:
    return sqrt(dot_2d(vector, vector))


@cy.profile(False)
@cy.ccall
def vector_magnitude_3d(vector: vec_xyz) -> cy.float:
    return sqrt(dot_3d(vector, vector))


# TODO: result uninitialized warning
@cy.profile(False)
@cy.ccall
def as_unit_vector_2d(vector: vec_xy) -> vec_xy:
    magnitude: cy.float = vector_magnitude_2d(vector)
    ux       : cy.float = vector[0] / magnitude
    uy       : cy.float = vector[1] / magnitude
    return (ux, uy)


@cy.profile(False)
@cy.ccall
def as_unit_vector_3d(vector: vec_xyz) -> vec_xyz:
    magnitude: cy.float = vector_magnitude_3d(vector)
    ux       : cy.float = vector[0] / magnitude
    uy       : cy.float = vector[1] / magnitude
    uz       : cy.float = vector[2] / magnitude
    return (ux, uy, uz)


# TODO: result uninitialized warning
@cy.profile(False)
@cy.ccall
def to_slope_plane(vector_2d: vec_xy, elevation_gradient: vec_xy) -> vec_xyz:
    return (
        vector_2d[0],
        vector_2d[1],
        dot_2d(vector_2d, elevation_gradient)
    )


# TODO: Replace numpy arrays with vec_xyz and vec_xy
@cy.profile(False)
@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def to_horizontal_plane(vector_3d: cy.double[:]) -> cy.double[:]:
    return vector_3d[0:2]


# TODO: Replace numpy arrays with vec_xyz and vec_xy
# TODO: Speed up conv.cartesian_to_azimuthal
@cy.profile(False)
@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def spread_direction_vector_to_angle(vector_3d: cy.double[:]) -> cy.float:
    vector_2d: cy.double[:] = to_horizontal_plane(vector_3d)
    x        : cy.double    = vector_2d[0]
    y        : cy.double    = vector_2d[1]
    az_coords: tuple        = conv.cartesian_to_azimuthal(x, y)
    azimuth  : cy.float     = az_coords[1]
    return azimuth


@cy.profile(False)
@cy.ccall
def get_slope_normal_vector(elevation_gradient: vec_xy) -> vec_xyz:
    slope_normal_vector: vec_xyz = (-elevation_gradient[0], -elevation_gradient[1], 1.0)
    return as_unit_vector_3d(slope_normal_vector)


# TODO: Eliminate need for numpy
# TODO: Speed up conv.cartesian_to_azimuthal and conv.opposite_direction
# TODO: Create a primitive math radians function
@cy.profile(False)
@cy.ccall
def rotate_on_sloped_plane(vector: np.ndarray, theta: cy.float, slope: cy.float, aspect: cy.float) -> vec_xyz:
    """
    Rotate a 3D vector <x,y,z> theta degrees clockwise on the plane defined by the slope and aspect.
    """
    # Calculate the slope normal vector from the slope and aspect
    elevation_gradient : vec_xy  = conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))
    slope_normal_vector: vec_xyz = get_slope_normal_vector(elevation_gradient)
    # Rotate theta degrees clockwise around the slope_normal_vector
    theta_rad: cy.float = radians(theta)
    return cos(theta_rad) * vector + np.cross(sin(theta_rad) * vector, slope_normal_vector)
# vector-utilities ends here
