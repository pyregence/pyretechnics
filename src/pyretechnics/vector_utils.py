# [[file:../../org/pyretechnics.org::vector-utilities][vector-utilities]]
# Fix cimport of numpy
import cython
if cython.compiled:
    from numpy import add, multiply, cross
    # from cython.cimports.numpy import add, multiply, cross
    from cython.cimports.pyretechnics.math import sqrt, sin, cos
    from cython.cimports.pyretechnics.cy_types import pyidx, vec_xy, vec_xyz
    from cython.cimports.pyretechnics.conversion import \
        opposite_direction, azimuthal_to_cartesian, cartesian_to_azimuthal, deg_to_rad
else:
    from numpy import add, multiply, cross
    from math import sqrt, sin, cos
    from pyretechnics.py_types import pyidx, vec_xy, vec_xyz
    from pyretechnics.conversion import \
        opposite_direction, azimuthal_to_cartesian, cartesian_to_azimuthal, deg_to_rad


import cython as cy

# TODO try inlining with @cy.inline (might require an inline directive in the .pxd file)


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
@cy.inline
@cy.exceptval(check=False)
def dot_3d(vector1: vec_xyz, vector2: vec_xyz) -> cy.float:
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]


@cy.profile(False)
@cy.ccall
def scale_2d(scalar: cy.float, vector: vec_xy) -> vec_xy:
    return (scalar * vector[0], scalar * vector[1])


@cy.profile(False)
@cy.ccall
def scale_3d(scalar: cy.float, vector: vec_xyz) -> vec_xyz:
    return (scalar * vector[0], scalar * vector[1], scalar * vector[2])

@cy.profile(False)
@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def add_2d(vector1: vec_xy, vector2: vec_xy) -> vec_xy:
    return (vector1[0] + vector2[0], vector1[1] + vector2[1])


@cy.profile(False)
@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def add_3d(vector1: vec_xyz, vector2: vec_xyz) -> vec_xyz:
    return (vector1[0] + vector2[0], vector1[1] + vector2[1], vector1[2] + vector2[2])


# TODO: result uninitialized warning
@cy.profile(False)
@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def vector_magnitude_2d(vector: vec_xy) -> cy.float:
    return sqrt(dot_2d(vector, vector))


@cy.profile(False)
@cy.ccall
@cy.inline
@cy.exceptval(check=False)
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
@cy.profile(False)
@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def spread_direction_vector_to_angle(vector_3d: vec_xyz) -> cy.float:
    x        : cy.double    = vector_3d[0]
    y        : cy.double    = vector_3d[1]
    az_coords: vec_xy       = cartesian_to_azimuthal(x, y)
    azimuth  : cy.float     = az_coords[1]
    return azimuth


@cy.profile(False)
@cy.ccall
def get_slope_normal_vector(elevation_gradient: vec_xy) -> vec_xyz:
    (dz_dx, dz_dy)               = elevation_gradient
    slope_normal_vector: vec_xyz = (-dz_dx, -dz_dy, 1.0)
    return as_unit_vector_3d(slope_normal_vector)


# TODO: Eliminate need for numpy
@cy.profile(False)
@cy.ccall
def rotate_on_sloped_plane(vector: cy.double[:], theta: cy.float, slope: cy.float, aspect: cy.float) -> cy.double[:]:
    """
    Rotate a 3D vector <x,y,z> theta degrees clockwise on the plane defined by the slope and aspect.
    """
    # Calculate the slope normal vector from the slope and aspect
    elevation_gradient : vec_xy  = azimuthal_to_cartesian(slope, opposite_direction(aspect))
    slope_normal_vector: vec_xyz = get_slope_normal_vector(elevation_gradient)
    # Rotate theta degrees clockwise around the slope_normal_vector
    theta_rad: cy.float = deg_to_rad(theta)
    return add(multiply(cos(theta_rad),
                        vector),
               cross(multiply(sin(theta_rad),
                              vector),
                    slope_normal_vector))
# vector-utilities ends here
