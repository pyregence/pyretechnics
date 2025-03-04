# [[file:../../org/pyretechnics.org::vector-utilities][vector-utilities]]
import cython
import cython as cy
if cython.compiled:
    from cython.cimports.libc.math import sqrt, sin, cos
    from cython.cimports.pyretechnics.cy_types import pyidx, vec_xy, vec_xyz
    import cython.cimports.pyretechnics.conversion as conv
else:
    from math import sqrt, sin, cos
    from pyretechnics.py_types import pyidx, vec_xy, vec_xyz
    import pyretechnics.conversion as conv


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def dot_2d(vector1: vec_xy, vector2: vec_xy) -> cy.float:
    return vector1[0] * vector2[0] + vector1[1] * vector2[1]


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def dot_3d(vector1: vec_xyz, vector2: vec_xyz) -> cy.float:
    return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def scale_2d(scalar: cy.float, vector: vec_xy) -> vec_xy:
    return (scalar * vector[0], scalar * vector[1])


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def scale_3d(scalar: cy.float, vector: vec_xyz) -> vec_xyz:
    return (scalar * vector[0], scalar * vector[1], scalar * vector[2])


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def add_2d(vector1: vec_xy, vector2: vec_xy) -> vec_xy:
    return (vector1[0] + vector2[0], vector1[1] + vector2[1])


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def add_3d(vector1: vec_xyz, vector2: vec_xyz) -> vec_xyz:
    return (vector1[0] + vector2[0], vector1[1] + vector2[1], vector1[2] + vector2[2])


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def vector_magnitude_2d(vector: vec_xy) -> cy.float:
    return sqrt(dot_2d(vector, vector))


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def vector_magnitude_3d(vector: vec_xyz) -> cy.float:
    return sqrt(dot_3d(vector, vector))


@cy.ccall
@cy.exceptval(check=False)
def as_unit_vector_2d(vector: vec_xy) -> vec_xy:
    magnitude: cy.float = vector_magnitude_2d(vector)
    if magnitude == 0.0:
        return vector
    else:
        ux: cy.float = vector[0] / magnitude
        uy: cy.float = vector[1] / magnitude
        return (ux, uy)


@cy.ccall
@cy.exceptval(check=False)
def as_unit_vector_3d(vector: vec_xyz) -> vec_xyz:
    magnitude: cy.float = vector_magnitude_3d(vector)
    if magnitude == 0.0:
        return vector
    else:
        ux: cy.float = vector[0] / magnitude
        uy: cy.float = vector[1] / magnitude
        uz: cy.float = vector[2] / magnitude
        return (ux, uy, uz)


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def to_slope_plane(vector_2d: vec_xy, elevation_gradient: vec_xy) -> vec_xyz:
    return (
        vector_2d[0],
        vector_2d[1],
        dot_2d(vector_2d, elevation_gradient)
    )


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def to_horizontal_plane(vector_3d: vec_xyz) -> vec_xy:
    return (vector_3d[0], vector_3d[1])


@cy.ccall
@cy.exceptval(check=False)
def spread_direction_vector_to_angle(vector_3d: vec_xyz) -> cy.float:
    x        : cy.float = vector_3d[0]
    y        : cy.float = vector_3d[1]
    az_coords: vec_xy   = conv.cartesian_to_azimuthal(x, y)
    azimuth  : cy.float = az_coords[1]
    return azimuth


@cy.ccall
@cy.exceptval(check=False)
def get_slope_normal_vector(elevation_gradient: vec_xy) -> vec_xyz:
    (dz_dx, dz_dy)               = elevation_gradient
    slope_normal_vector: vec_xyz = (-dz_dx, -dz_dy, 1.0)
    return as_unit_vector_3d(slope_normal_vector)


@cy.ccall
@cy.exceptval(check=False)
def cross_3d(vector1: vec_xyz, vector2: vec_xyz) -> vec_xyz:
    (a, b, c) = vector1
    (d, e, f) = vector2
    return (
        b * f - e * c,
        -a * f + d * c,
        a * e - d * b,
    )


@cy.ccall
@cy.exceptval(check=False)
def rotate_on_sloped_plane(vector: vec_xyz, theta: cy.float, slope: cy.float, aspect: cy.float) -> vec_xyz:
    """
    Rotate a 3D vector <x,y,z> theta degrees clockwise on the plane defined by the slope and aspect.
    """
    # Calculate the slope normal vector from the slope and aspect
    elevation_gradient : vec_xy  = conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))
    slope_normal_vector: vec_xyz = get_slope_normal_vector(elevation_gradient)
    # Calculate sine and cosine of theta
    theta_rad: cy.float = conv.deg_to_rad(theta)
    cos_theta: cy.float = cos(theta_rad)
    sin_theta: cy.float = sin(theta_rad)
    # Rotate theta degrees clockwise around the slope_normal_vector
    vector_i: vec_xyz = (
        cos_theta * vector[0],
        cos_theta * vector[1],
        cos_theta * vector[2],
    )
    vector_j: vec_xyz = (
        sin_theta * vector[0],
        sin_theta * vector[1],
        sin_theta * vector[2],
    )
    vector_k: vec_xyz = cross_3d(vector_j, slope_normal_vector)
    return (
        vector_i[0] + vector_k[0],
        vector_i[1] + vector_k[1],
        vector_i[2] + vector_k[2],
    )
# vector-utilities ends here
