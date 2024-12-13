# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients-approx][phi-field-spatial-gradients-approx]]
import cython
if cython.compiled:
    from cython.cimports.cpython.mem import PyMem_Malloc, PyMem_Realloc, PyMem_Free
    from cython.cimports.pyretechnics.math import sqrt, atan
    from cython.cimports.pyretechnics.cy_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, fcatarr, fclaarr, FuelModel, FireBehaviorMax, SpreadBehavior
    import cython.cimports.pyretechnics.crown_fire as cf
    from cython.cimports.pyretechnics.conversion import \
        rad_to_deg, opposite_direction, azimuthal_to_cartesian, wind_speed_10m_to_wind_speed_20ft, \
        Btu_lb_to_kJ_kg, km_hr_to_m_min, m_to_ft
    import cython.cimports.pyretechnics.narrow_band_tracking as nbt
    import cython.cimports.pyretechnics.surface_fire1 as sf
    from cython.cimports.pyretechnics.space_time_cube import ISpaceTimeCube
    from cython.cimports.pyretechnics.random import BufferedRandGen
    import cython.cimports.pyretechnics.spot_fire as spot
    from cython.cimports.pyretechnics.vector_utils import \
        vector_magnitude_2d, vector_magnitude_3d, as_unit_vector_2d, as_unit_vector_3d, dot_2d, dot_3d, scale_2d, scale_3d, \
        get_slope_normal_vector, to_slope_plane, spread_direction_vector_to_angle
else:
    from math import sqrt, atan
    from pyretechnics.py_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, fcatarr, fclaarr, FuelModel, FireBehaviorMax, SpreadBehavior
    from pyretechnics.conversion import \
        rad_to_deg, opposite_direction, azimuthal_to_cartesian, wind_speed_10m_to_wind_speed_20ft, \
        Btu_lb_to_kJ_kg, km_hr_to_m_min, m_to_ft
    import pyretechnics.crown_fire as cf
    from pyretechnics.space_time_cube import ISpaceTimeCube
    import pyretechnics.spot_fire as spot
    import pyretechnics.narrow_band_tracking as nbt
    from pyretechnics.random import BufferedRandGen
    import pyretechnics.surface_fire1 as sf
    from pyretechnics.vector_utils import \
        vector_magnitude_2d, vector_magnitude_3d, as_unit_vector_2d, as_unit_vector_3d, dot_2d, dot_3d, scale_2d, scale_3d, \
        get_slope_normal_vector, to_slope_plane, spread_direction_vector_to_angle


import cython as cy
# TODO: cimport all of the modules below
import numpy as np
import pyretechnics.fuel_models as fm
import pyretechnics.surface_fire1 as sf0
import sortedcontainers as sortc

PI = cy.declare(cy.double, 3.14159265358979323846)


@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_dphi_dx_approx(phi: cy.float[:,:], dx: cy.float, x: pyidx, y: pyidx, cols: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given the cell width dx.
    """
    east_x: pyidx = x + 1
    west_x: pyidx = x - 1
    if east_x < cols:
        if west_x >= 0:
            return (phi[y][east_x] - phi[y][west_x]) / (2.0 * dx)
        else:
            return (phi[y][east_x] - phi[y][x]) / dx
    else:
        if west_x >= 0:
            return (phi[y][x] - phi[y][west_x]) / dx
        else:
            return 0.0


@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_dphi_dy_approx(phi: cy.float[:,:], dy: cy.float, x: pyidx, y: pyidx, rows: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given the cell height dy.
    """
    north_y: pyidx = y + 1
    south_y: pyidx = y - 1
    if north_y < rows:
        if south_y >= 0:
            return (phi[north_y][x] - phi[south_y][x]) / (2.0 * dy)
        else:
            return (phi[north_y][x] - phi[y][x]) / dy
    else:
        if south_y >= 0:
            return (phi[y][x] - phi[south_y][x]) / dy
        else:
            return 0.0


# TODO: Pass rows and cols
# TODO: Handle exception values from child functions
@cy.profile(False)
@cy.ccall
def calc_phi_gradient_approx(phi: cy.float[:,:], dx: cy.float, dy: cy.float, x: pyidx, y: pyidx) -> vec_xy:
    """
    Calculate the spatial gradient of the phi raster at grid cell (x,y)
    given the cell width dx and the cell height dy.
    """
    rows   : pyidx    = phi.shape[0]
    cols   : pyidx    = phi.shape[1]
    dphi_dx: cy.float = calc_dphi_dx_approx(phi, dx, x, y, cols)
    dphi_dy: cy.float = calc_dphi_dy_approx(phi, dy, x, y, rows)
    return (dphi_dx, dphi_dy)
# phi-field-spatial-gradients-approx ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector][phi-field-normal-vector]]
# TODO: Remove unused function
@cy.profile(False)
@cy.ccall
def calc_phi_normal_vector(phi: cy.float[:,:], dx: cy.float, dy: cy.float, x: pyidx, y: pyidx) -> vec_xy:
    """
    Calculate the phi field normal vector in the x and y dimensions.

    - n_x: eastward component of the unit normal vector
    - n_y: northward component of the unit normal vector
    """
    phi_gradient: vec_xy = calc_phi_gradient_approx(phi, dx, dy, x, y)
    if phi_gradient[0] == 0.0 and phi_gradient[1] == 0.0:
        return phi_gradient # (n_x, n_y)
    else:
        return as_unit_vector_2d(phi_gradient) # (n_x, n_y)
# phi-field-normal-vector ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector-angle][phi-field-normal-vector-angle]]
# TODO: Remove unused function
@cy.profile(False)
@cy.ccall
@cy.cdivision(True)
def calc_phi_normal_azimuth(phi_normal_vector: vec_xy) -> cy.float:
    """
    Calculate the angle (measured in degrees clockwise from North)
    to which the phi field's normal vector points.
    """
    (n_x, n_y) = phi_normal_vector
    angle: cy.float
    if n_x > 0.0:
        if n_y >= 0.0:
            angle = 0.5 * PI - atan(n_y / n_x)
        elif n_y < 0.0:
            angle = 0.5 * PI + atan(abs(n_y) / n_x)
    elif n_x < 0.0:
        if n_y >= 0.0:
            angle = 1.5 * PI + atan(n_y / abs(n_x))
        elif n_y < 0.0:
            angle = 1.5 * PI - atan(n_y / n_x)
    else:
        if n_y >= 0.0:
            angle = 0.0
        elif n_y < 0.0:
            angle = PI
    return rad_to_deg(angle)
# phi-field-normal-vector-angle ends here
# [[file:../../org/pyretechnics.org::superbee-flux-limiter][superbee-flux-limiter]]
@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
def calc_superbee_flux_limiter(dphi_up: cy.float, dphi_loc: cy.float) -> cy.float:
    """
    TODO: Add docstring
    """
    if dphi_loc == 0.0:
        return 0.0
    else:
        r: cy.float = dphi_up / dphi_loc
        return max(0.0,
                   min(2.0 * r, 1.0),
                   min(r, 2.0))
# superbee-flux-limiter ends here
# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients][phi-field-spatial-gradients]]
@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
def calc_dphi_dx(phi: cy.float[:,:], u_x: cy.float, dx: cy.float, x: pyidx, y: pyidx, cols: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_x  :: m/min
    - dx   :: meters
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - cols :: integer number of columns in the phi matrix
    """
    phi_east: cy.float = calc_phi_east(phi, u_x, x, y, cols)
    phi_west: cy.float = calc_phi_west(phi, u_x, x, y, cols)
    return (phi_east - phi_west) / dx


@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
def calc_dphi_dy(phi: cy.float[:,:], u_y: cy.float, dy: cy.float, x: pyidx, y: pyidx, rows: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_y  :: m/min
    - dy   :: meters
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - rows :: integer number of rows in the phi matrix
    """
    phi_north: cy.float = calc_phi_north(phi, u_y, x, y, rows)
    phi_south: cy.float = calc_phi_south(phi, u_y, x, y, rows)
    return (phi_north - phi_south) / dy

# TODO: Pass rows and cols
# TODO: Handle exception values from child functions
@cy.profile(False)
@cy.ccall
def calc_phi_gradient(phi: cy.float[:,:], u_x: cy.float, u_y: cy.float, dx: cy.float, dy: cy.float,
                      x: pyidx, y: pyidx) -> vec_xy:
    """
    Calculate the spatial gradient of the phi raster at grid cell (x,y) given:
    - phi :: 2D float array of values in [-1,1]
    - u_x :: m/min
    - u_y :: m/min
    - dx  :: meters
    - dy  :: meters
    - x   :: integer column index in phi
    - y   :: integer row index in phi
    """
    rows   : pyidx    = phi.shape[0]
    cols   : pyidx    = phi.shape[1]
    dphi_dx: cy.float = calc_dphi_dx(phi, u_x, dx, x, y, cols)
    dphi_dy: cy.float = calc_dphi_dy(phi, u_y, dy, x, y, rows)
    return (dphi_dx, dphi_dy)
# phi-field-spatial-gradients ends here
# [[file:../../org/pyretechnics.org::phi-east][phi-east]]
@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_east(phi: cy.float[:,:], u_x: cy.float, x: pyidx, y: pyidx, cols: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_x  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - cols :: integer number of columns in the phi matrix
    """
    very_east_x: pyidx = min(x+2, cols-1)
    east_x     : pyidx = min(x+1, cols-1)
    west_x     : pyidx = max(x-1, 0)

    dphi_loc: cy.float = phi[y][east_x] - phi[y][x] # FIXME bad array lookup notation I think?
    if u_x >= 0.0:
        dphi_up: cy.float = phi[y][x] - phi[y][west_x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        dphi_up: cy.float = phi[y][very_east_x] - phi[y][east_x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][east_x] - 0.5 * B * dphi_loc
# phi-east ends here
# [[file:../../org/pyretechnics.org::phi-west][phi-west]]
@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_west(phi: cy.float[:,:], u_x: cy.float, x: pyidx, y: pyidx, cols: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the -x (east->west)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_x  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - cols :: integer number of columns in the phi matrix
    """
    east_x     : pyidx = min(x+1, cols-1)
    west_x     : pyidx = max(x-1, 0)
    very_west_x: pyidx = max(x-2, 0)

    dphi_loc: cy.float = phi[y][west_x] - phi[y][x]
    if u_x >= 0.0:
        dphi_up: cy.float = phi[y][very_west_x] - phi[y][west_x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][west_x] - 0.5 * B * dphi_loc
    else:
        dphi_up: cy.float = phi[y][x] - phi[y][east_x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
# phi-west ends here
# [[file:../../org/pyretechnics.org::phi-north][phi-north]]
@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_north(phi: cy.float[:,:], u_y: cy.float, x: pyidx, y: pyidx, rows: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_y  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - rows :: integer number of rows in the phi matrix
    """
    very_north_y: pyidx = min(y+2, rows-1)
    north_y     : pyidx = min(y+1, rows-1)
    south_y     : pyidx = max(y-1, 0)

    dphi_loc: cy.float = phi[north_y][x] - phi[y][x]
    if u_y >= 0.0:
        dphi_up: cy.float = phi[y][x] - phi[south_y][x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
    else:
        dphi_up: cy.float = phi[very_north_y][x] - phi[north_y][x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[north_y][x] - 0.5 * B * dphi_loc
# phi-north ends here
# [[file:../../org/pyretechnics.org::phi-south][phi-south]]
@cy.profile(False)
@cy.cfunc
@cy.exceptval(65504.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_south(phi: cy.float[:,:], u_y: cy.float, x: pyidx, y: pyidx, rows: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the -y (north->south)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_y  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - rows :: integer number of rows in the phi matrix
    """
    north_y     : pyidx = min(y+1, rows-1)
    south_y     : pyidx = max(y-1, 0)
    very_south_y: pyidx = max(y-2, 0)

    dphi_loc: cy.float = phi[south_y][x] - phi[y][x]
    if u_y >= 0.0:
        dphi_up: cy.float = phi[very_south_y][x] - phi[south_y][x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[south_y][x] - 0.5 * B * dphi_loc
    else:
        dphi_up: cy.float = phi[y][x] - phi[north_y][x]
        B      : cy.float = calc_superbee_flux_limiter(dphi_up, dphi_loc)
        return phi[y][x] + 0.5 * B * dphi_loc
# phi-south ends here
# [[file:../../org/pyretechnics.org::calc-fireline-normal-behavior][calc-fireline-normal-behavior]]
# TODO: Move this to pyretechnics.vector_utils and use throughout the literate program
@cy.profile(True)
@cy.ccall
def calc_elevation_gradient(slope: cy.float, aspect: cy.float) -> vec_xy:
    """
    Returns the elevation gradient (dz_dx: rise/run, dz_dy: rise/run) given:
    - slope  :: rise/run
    - aspect :: degrees clockwise from North
    """
    return azimuthal_to_cartesian(slope, opposite_direction(aspect))

# NOTE it would be better to use a cython enum here, but that's not supported in pure python syntax...
fire_type_unburned      = cy.declare(cy.int, 0)
fire_type_surface       = cy.declare(cy.int, 1)
fire_type_crown_passive = cy.declare(cy.int, 2)
fire_type_crown_active  = cy.declare(cy.int, 3)

def stringify_fire_type(fire_type: integer) -> string: # TODO only here for backwards compatibility
    return ("unburned" if fire_type == fire_type_unburned
            else "surface" if fire_type == fire_type_surface
            else "passive_crown" if fire_type == fire_type_crown_passive
            else "active_crown")


@cy.profile(True)
@cy.ccall
def calc_phi_gradient_on_slope(phi_gradient_xy: vec_xy, elevation_gradient: vec_xy) -> vec_xyz:
    """
    Return the gradient of phi projected onto the slope-tangential plane as a 3D (x,y,z) vector (in phi/m) given:
    - phi_gradient_xy    :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - elevation_gradient :: (dz_dx: rise/run, dz_dy: rise/run)
    """
    (dphi_dx, dphi_dy)        = phi_gradient_xy
    phi_gradient_xyz: vec_xyz = (dphi_dx, dphi_dy, 0.0)
    if vector_magnitude_2d(elevation_gradient) == 0.0:
        return phi_gradient_xyz
    else:
        slope_normal_vector: vec_xyz  = get_slope_normal_vector(elevation_gradient) # (x,y,z) unit vector
        phi_slope_agreement: cy.float = dot_3d(phi_gradient_xyz, slope_normal_vector)
        dphi_dx_on_slope   : cy.float = phi_gradient_xyz[0] - phi_slope_agreement * slope_normal_vector[0]
        dphi_dy_on_slope   : cy.float = phi_gradient_xyz[1] - phi_slope_agreement * slope_normal_vector[1]
        dphi_dz_on_slope   : cy.float = phi_gradient_xyz[2] - phi_slope_agreement * slope_normal_vector[2]
        return (dphi_dx_on_slope, dphi_dy_on_slope, dphi_dz_on_slope)


# FIXME: Do I switch to cruz_passive_crown_fire_spread_rate() if the normal_spread_rate < critical_spread_rate?
#        Did I do this correctly in calc_crown_fire_behavior_in_direction?
@cy.profile(True)
@cy.ccall
def calc_fireline_normal_behavior(fire_behavior_max: FireBehaviorMax, phi_gradient: vec_xyz) -> SpreadBehavior:
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

    phi_magnitude = vector_magnitude_3d(phi_gradient) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burning
    #================================================================================================

    if (phi_magnitude == 0.0 or fire_behavior_max.max_spread_rate == 0.0):
        # This location is not on the fire perimeter and/or is not burning

        #================================================================================================
        # Set the spread direction to the phi gradient direction, max spread direction, upslope, or North
        #================================================================================================

        spread_direction = (np.asarray(phi_gradient) / phi_magnitude if phi_magnitude > 0.0
                            else fire_behavior_max.max_spread_direction)
        

        #============================================================================================
        # Return zero surface/crown fire behavior
        #============================================================================================

        sd: vec_xyz = (spread_direction[0], spread_direction[1], spread_direction[2])
        return SpreadBehavior(
            dphi_dt            = 0.0,
            fire_type          = fire_type_unburned,
            spread_rate        = 0.0,
            spread_direction   = sd,
            fireline_intensity = 0.0,
            flame_length       = 0.0
        )

    else:
        # This location is on the fire perimeter and is burning

        #============================================================================================
        # Unpack the fire_behavior_max dictionary
        #============================================================================================

        heading_fire_type: cy.int = fire_behavior_max.max_fire_type # FIXME .get("max_fire_type", "surface")
        if heading_fire_type < 0:
            heading_fire_type = fire_type_surface
        heading_spread_rate: cy.float        = fire_behavior_max.max_spread_rate               # m/min
        heading_spread_direction: vec_xyz    = fire_behavior_max.max_spread_direction          # (x,y,z) unit vector
        heading_spread_vector: vec_xyz       = scale_3d(heading_spread_rate, heading_spread_direction)     # (x,y,z) m/min vector
        heading_fireline_intensity: cy.float = fire_behavior_max.max_fireline_intensity        # kW/m
        length_to_width_ratio: cy.float      = fire_behavior_max.length_to_width_ratio         # unitless
        eccentricity: cy.float               = fire_behavior_max.eccentricity                  # unitless
        critical_spread_rate: cy.float       = fire_behavior_max.critical_spread_rate          # m/min
        if critical_spread_rate < 0.0:
            critical_spread_rate = 0.0

        #============================================================================================
        # Calculate the backing and flanking fire spread rates
        #============================================================================================

        backing_adjustment: cy.float   = (1.0 - eccentricity) / (1.0 + eccentricity)                                 # unitless
        backing_spread_rate: cy.float  = heading_spread_rate * backing_adjustment                                    # m/min
        flanking_spread_rate: cy.float = (heading_spread_rate + backing_spread_rate) / (2.0 * length_to_width_ratio) # m/min

        #============================================================================================
        # Calculate dphi/dt
        #============================================================================================

        A: cy.float = (heading_spread_rate - backing_spread_rate) / (2 * heading_spread_rate) # unitless

        B: cy.float = dot_3d(heading_spread_vector, phi_gradient)                             # phi/min
        C: cy.float = flanking_spread_rate / heading_spread_rate                              # unitless
        D: cy.float = (heading_spread_rate * phi_magnitude) ** 2.0                            # (phi/min)^2
        E: cy.float = (length_to_width_ratio ** 2.0 - 1.0) * (B ** 2.0)                       # (phi/min)^2
        dphi_dt: cy.float = -(A * B + C * sqrt(D + E))                                              # phi/min

        #============================================================================================
        # Calculate fire behavior normal to the fire perimeter
        #============================================================================================

        normal_spread_rate: cy.float        = -dphi_dt / phi_magnitude                        # m/min
        normal_direction: vec_xyz           = as_unit_vector_3d(phi_gradient)                 # (x,y,z) unit vector
        normal_adjustment: cy.float         = normal_spread_rate / heading_spread_rate        # unitless
        normal_fireline_intensity: cy.float = heading_fireline_intensity * normal_adjustment  # kW/m
        normal_flame_length: cy.float = sf.calc_flame_length(normal_fireline_intensity) # m
        normal_fire_type: cy.int = (
            fire_type_surface if heading_fire_type == fire_type_surface
            else fire_type_crown_active if normal_spread_rate > critical_spread_rate
            else fire_type_crown_passive)

        #========================================================================================
        # Return the surface/crown fire behavior normal to the fire perimeter
        #========================================================================================

        sd: vec_xyz = (normal_direction[0], normal_direction[1], normal_direction[2])
        return SpreadBehavior(
            dphi_dt            = dphi_dt,                   # phi/min
            fire_type          = normal_fire_type,          # surface, passive_crown, or active_crown
            spread_rate        = normal_spread_rate,        # m/min
            spread_direction   = sd,          # (x,y,z) unit vector
            fireline_intensity = normal_fireline_intensity, # kW/m
            flame_length       = normal_flame_length        # m
        )
# calc-fireline-normal-behavior ends here
# [[file:../../org/pyretechnics.org::burn-cell-toward-phi-gradient][burn-cell-toward-phi-gradient]]
# TODO: Create a version of this function that runs efficiently over a space_time_region

@cy.cfunc
def fclaarr_from_list(l: object) -> fclaarr:
    l0: cy.float = l[0]
    l1: cy.float = l[1]
    l2: cy.float = l[2]
    l3: cy.float = l[3]
    l4: cy.float = l[4]
    l5: cy.float = l[5]
    ret: fclaarr = (l0, l1, l2, l3, l4, l5)
    return ret

@cy.cfunc
@cy.cdivision(True)
def compute_f_A_sigma(A: cy.float, sigma_ij: cy.float) -> cy.float:
    if sigma_ij > 0:
        return exp(A / sigma_ij)
    else:
        return 0


@cy.cfunc
@cy.profile(False)
@cy.exceptval(check=False)
def firemod_size_class(sigma_i: cy.float) -> cy.float:
    s: cy.float = sigma_i    
    return (
        1 if (s >= 1200.0)
        else 2 if (s >= 192.0)
        else 3 if (s >= 96.0)
        else 4 if (s >= 48.0)
        else 5 if (s >= 16.0)
        else 6
    )

@cy.profile(False)
@cy.ccall
def fm_struct(fm: dict) -> FuelModel:
    sigma: fclaarr = fclaarr_from_list(fm["sigma"])
    return FuelModel(
        number   = fm["number"],
        delta    = fm["delta"],
        M_x      = fclaarr_from_list(fm["M_x"]),
        M_f      = fclaarr_from_list(fm["M_x"]),
        w_o      = fclaarr_from_list(fm["w_o"]),
        sigma    = sigma,
        h        = fclaarr_from_list(fm["h"]),
        rho_p    = fclaarr_from_list(fm["rho_p"]),
        S_T      = fclaarr_from_list(fm["S_T"]),
        S_e      = fclaarr_from_list(fm["S_e"]),
        dynamic  = fm["dynamic"],
        burnable = fm["burnable"],
        exp_A_sigma = (
            compute_f_A_sigma(-138.0, sigma[0]),
            compute_f_A_sigma(-138.0, sigma[1]),
            compute_f_A_sigma(-138.0, sigma[2]),
            compute_f_A_sigma(-138.0, sigma[3]),
            compute_f_A_sigma(-500.0, sigma[4]),
            compute_f_A_sigma(-500.0, sigma[5])
        ),
        firemod_size_classes = (
            firemod_size_class(sigma[0]),
            firemod_size_class(sigma[1]),
            firemod_size_class(sigma[2]),
            firemod_size_class(sigma[3]),
            firemod_size_class(sigma[4]),
            firemod_size_class(sigma[5]),
        ),
        f_ij    = (0.,0.,0.,0.,0.,0.),
        f_i = (0.,0.),
        g_ij = (0.,0.,0.,0.,0.,0.)
    )

if cython.compiled:
    from cython.cimports.pyretechnics.math import exp
else:
    from math import exp

@cy.cfunc
@cy.profile(False)
#@cy.inline
def add_dynamic_fuel_loading(fuel_model: FuelModel, M_f: fclaarr) -> FuelModel:
    """
    Updates M_f and w_o.
    """
    if fuel_model.dynamic:
        # dynamic fuel model
        w_o: fclaarr              = fuel_model.w_o
        live_herbaceous_load: cy.float      = w_o[4]
        live_herbaceous_moisture: cy.float  = M_f[4]
        fraction_green: cy.float            = max(0.0, min(1.0, (live_herbaceous_moisture / 0.9) - 0.3333333333333333))
        fraction_cured: cy.float            = 1.0 - fraction_green
        dynamic_fuel_model        = fuel_model
        M_f1: fclaarr = M_f
        M_f1[3] = M_f1[0] # set dead_herbaceous to dead_1hr
        dynamic_fuel_model.M_f = M_f1
        w_o1: fclaarr = w_o
        w_o1[3] = live_herbaceous_load * fraction_cured # dead_herbaceous
        w_o1[4] = live_herbaceous_load * fraction_green # live_herbaceous
        dynamic_fuel_model.w_o = w_o1
        return dynamic_fuel_model
    else:
        # static fuel model
        static_fuel_model = fuel_model
        static_fuel_model.M_f = M_f
        return static_fuel_model


@cy.cfunc
@cy.profile(False)
@cy.inline
@cy.exceptval(check=False)
@cy.cdivision(True)
def _compute_f_ij(A_ij: cy.float, A: cy.float) -> cy.float:
    return (A_ij / A) if A > 0.0 else 0.0




@cy.cfunc
@cy.profile(False)
@cy.exceptval(check=False)
def gij_i(firemod_size_classes: fclaarr, f_ij: fclaarr, firemod_size_class_i: cy.float, is_dead: cy.bint) -> cy.float:
    """
    Sums the f_ij of the same category (dead/live) as i, and having the same firemod_size_class.
    """
    c: cy.float = firemod_size_class_i
    # NOTE there may be repetitions in firemod_size_classes, this is why this expression is not trivially equal to f_ij[i]:
    return (( # TODO OPTIM pre-compute this conditional branching (it's fully determined by sigma). Might be represented efficienty in bit flags.
        (f_ij[0] if (c == firemod_size_classes[0]) else 0.0)
        + (f_ij[1] if (c == firemod_size_classes[1]) else 0.0)
            + (f_ij[2] if (c == firemod_size_classes[2]) else 0.0)
            + (f_ij[3] if (c == firemod_size_classes[3]) else 0.0))
        if is_dead else
        ((f_ij[4] if (c == firemod_size_classes[4]) else 0.0)
            + (f_ij[5] if (c == firemod_size_classes[5]) else 0.0)))


@cy.cfunc
@cy.cdivision(True)
#@cy.inline
@cy.profile(False)
def add_weighting_factors(fuel_model: FuelModel) -> FuelModel:
    """
    Assigns f_ij, f_i and g_ij.
    
    """
    w_o: fclaarr                         = fuel_model.w_o
    sigma: fclaarr                       = fuel_model.sigma
    rho_p: fclaarr                       = fuel_model.rho_p
    A_ij: fclaarr = ( # TODO OPTIM pre-compute sigma/rho_p
        (sigma[0] * w_o[0]) / rho_p[0],
        (sigma[1] * w_o[1]) / rho_p[1],
        (sigma[2] * w_o[2]) / rho_p[2],
        (sigma[3] * w_o[3]) / rho_p[3],
        (sigma[4] * w_o[4]) / rho_p[4],
        (sigma[5] * w_o[5]) / rho_p[5]
    )
    A_i: fcatarr = (
        A_ij[0] + A_ij[1] + A_ij[2] + A_ij[3],
        A_ij[4] + A_ij[5]
    )
    A_0, A_1 = A_i
    f_ij: fclaarr = (0, 0, 0, 0, 0, 0)
    if A_0 > 0:
        A0_inv: cy.float = 1. / A_0
        f_ij[0] = A_ij[0] * A0_inv
        f_ij[1] = A_ij[1] * A0_inv
        f_ij[2] = A_ij[2] * A0_inv
        f_ij[3] = A_ij[3] * A0_inv
    if A_1 > 0:
        A1_inv: cy.float = 1. / A_1
        f_ij[4] = A_ij[4] * A1_inv
        f_ij[5] = A_ij[5] * A1_inv
    A_T: cy.float = A_0 + A_1
    f_i: fcatarr = (0, 0)
    if A_T > 0.0:
        A_T_inv: cy.float = 1. / A_T
        f_i[0] = A_0 * A_T_inv
        f_i[1] = A_1 * A_T_inv
    firemod_size_classes: fclaarr = fuel_model.firemod_size_classes
    g_ij: fclaarr = (
        gij_i(firemod_size_classes, f_ij, firemod_size_classes[0], True),
        gij_i(firemod_size_classes, f_ij, firemod_size_classes[1], True),
        gij_i(firemod_size_classes, f_ij, firemod_size_classes[2], True),
        gij_i(firemod_size_classes, f_ij, firemod_size_classes[3], True),
        gij_i(firemod_size_classes, f_ij, firemod_size_classes[4], False),
        gij_i(firemod_size_classes, f_ij, firemod_size_classes[5], False)
    )
    weighted_fuel_model      = fuel_model
    weighted_fuel_model.f_ij = f_ij
    weighted_fuel_model.f_i  = f_i
    weighted_fuel_model.g_ij = g_ij
    return weighted_fuel_model

@cy.profile(False)
@cy.cdivision(True)
@cy.cfunc
#@cy.inline
@cy.exceptval(check=False)
def _lf(sigma_ij: cy.float, i: cy.int, w_o_i: cy.float) -> cy.float:
    A: cy.float = -138.0 if (i < 4) else -500.0
    lfi: cy.float
    if (sigma_ij > 0.0):
        lfi = w_o_i * exp(A / sigma_ij) # TODO OPTIM pre-compute the exponential.
    else:
        lfi = 0.0
    return lfi


@cy.cfunc
@cy.cdivision(True)
#@cy.inline
@cy.profile(False)
def add_live_moisture_of_extinction(fuel_model: FuelModel) -> FuelModel:
    """
    Equation 88 from Rothermel 1972 adjusted by Albini 1976 Appendix III.

    Updates M_x.
    """
    w_o: fclaarr                       = fuel_model.w_o
    sigma: fclaarr                     = fuel_model.sigma
    M_f: fclaarr                       = fuel_model.M_f
    M_x: fclaarr                       = fuel_model.M_x
    exp_A_sigma: fclaarr               = fuel_model.exp_A_sigma
    loading_factors: fclaarr = (
        w_o[0] * exp_A_sigma[0],
        w_o[1] * exp_A_sigma[1],
        w_o[2] * exp_A_sigma[2],
        w_o[3] * exp_A_sigma[3],
        w_o[4] * exp_A_sigma[4],
        w_o[5] * exp_A_sigma[5]
    )
    lf: fclaarr = loading_factors
    dead_loading_factor: cy.float = lf[0] + lf[1] + lf[2] + lf[3]
    live_loading_factor: cy.float = lf[4] + lf[5]
    dead_moisture_factor: cy.float = (
        M_f[0] * loading_factors[0] +
        M_f[1] * loading_factors[1] +
        M_f[2] * loading_factors[2] +
        M_f[3] * loading_factors[3]
    )
    dead_fuel_moisture: cy.float        = (dead_moisture_factor / dead_loading_factor) if (dead_loading_factor > 0.0) else 0.0
    M_x_dead: cy.float = M_x[0]
    M_x_live: cy.float
    if (live_loading_factor > 0.0):
        dead_to_live_ratio: cy.float = (dead_loading_factor / live_loading_factor)
        M_x_live = max(
            M_x_dead,
            (2.9 * dead_to_live_ratio * (1.0 - (dead_fuel_moisture / M_x_dead))) - 0.226)
    else:
        M_x_live = M_x_dead
    moisturized_fuel_model    = fuel_model
    M_x1: fclaarr = (
        M_x[0],
        M_x[1],
        M_x[2],
        M_x[3],
        M_x_live,
        M_x_live,
    )
    moisturized_fuel_model.M_x = M_x1
    return moisturized_fuel_model

@cy.cfunc
def moisturize(fuel_model: FuelModel, fuel_moisture: fclaarr) -> FuelModel:
    """
    Updates w_o and M_f, and assigns M_x, f_ij, f_i and g_ij.
    """
    dynamic_fuel_model     = add_dynamic_fuel_loading(fuel_model, fuel_moisture)
    weighted_fuel_model    = add_weighting_factors(dynamic_fuel_model)
    moisturized_fuel_model = add_live_moisture_of_extinction(weighted_fuel_model)
    return moisturized_fuel_model
# moisturize ends here

def fm_structs() -> dict[int, FuelModel]:
    return {f["number"]: fm_struct(f) for f in fm.fuel_model_table.values()}


fuel_model_structs = fm_structs()



@cy.cclass
class SpreadInputs:
    """
    A fast-access data structure for reading inputs in performance-critical code.
    """
    rows: pyidx
    cols: pyidx
    band_duration: cy.float # seconds
    spatial_resolution: vec_xy # meters
    

    slope: ISpaceTimeCube
    aspect: ISpaceTimeCube
    fuel_model: ISpaceTimeCube
    canopy_cover: ISpaceTimeCube
    canopy_height: ISpaceTimeCube
    canopy_base_height: ISpaceTimeCube
    canopy_bulk_density: ISpaceTimeCube
    wind_speed_10m: ISpaceTimeCube
    upwind_direction: ISpaceTimeCube
    fuel_moisture_dead_1hr: ISpaceTimeCube
    fuel_moisture_dead_10hr: ISpaceTimeCube
    fuel_moisture_dead_100hr: ISpaceTimeCube
    fuel_moisture_live_herbaceous: ISpaceTimeCube
    fuel_moisture_live_woody: ISpaceTimeCube
    foliar_moisture: ISpaceTimeCube
    temperature: ISpaceTimeCube
    fuel_spread_adjustment: ISpaceTimeCube
    weather_spread_adjustment: ISpaceTimeCube
    
    fuel_models_arr: cy.pointer(FuelModel)

    def __cinit__(self, cube_resolution,
                 slope: ISpaceTimeCube,
                 aspect: ISpaceTimeCube,
                 fuel_model: ISpaceTimeCube,
                 canopy_cover: ISpaceTimeCube,
                 canopy_height: ISpaceTimeCube,
                 canopy_base_height: ISpaceTimeCube,
                 canopy_bulk_density: ISpaceTimeCube,
                 wind_speed_10m: ISpaceTimeCube,
                 upwind_direction: ISpaceTimeCube,
                 fuel_moisture_dead_1hr: ISpaceTimeCube,
                 fuel_moisture_dead_10hr: ISpaceTimeCube,
                 fuel_moisture_dead_100hr: ISpaceTimeCube,
                 fuel_moisture_live_herbaceous: ISpaceTimeCube,
                 fuel_moisture_live_woody: ISpaceTimeCube,
                 foliar_moisture: ISpaceTimeCube,
                 temperature: ISpaceTimeCube,
                 fuel_spread_adjustment: ISpaceTimeCube,
                 weather_spread_adjustment: ISpaceTimeCube
                 ):
        self.band_duration   = cube_resolution[0]
        (_bands, rows, cols) = slope.shape
        self.rows = rows
        self.cols = cols
        self.spatial_resolution = (rows, cols)

        self.slope = slope
        self.aspect = aspect
        self.fuel_model = fuel_model
        self.canopy_cover = canopy_cover
        self.canopy_height = canopy_height
        self.canopy_base_height = canopy_base_height
        self.canopy_bulk_density = canopy_bulk_density
        self.wind_speed_10m = wind_speed_10m
        self.upwind_direction = upwind_direction
        self.fuel_moisture_dead_1hr = fuel_moisture_dead_1hr
        self.fuel_moisture_dead_10hr = fuel_moisture_dead_10hr
        self.fuel_moisture_dead_100hr = fuel_moisture_dead_100hr
        self.fuel_moisture_live_herbaceous = fuel_moisture_live_herbaceous
        self.fuel_moisture_live_woody = fuel_moisture_live_woody
        self.foliar_moisture = foliar_moisture
        self.temperature = temperature
        self.fuel_spread_adjustment = fuel_spread_adjustment
        self.weather_spread_adjustment = weather_spread_adjustment
        self._init_fuel_models()


    @cy.ccall
    def _init_fuel_models(self):
        # fm_arr = cvarray(shape=(300,), itemsize=cython.sizeof(FuelModel), format="i")
        # fuel_models_arr = cython.declare(FuelModel[:], fm_arr)
        self.fuel_models_arr = cy.cast(
            cy.pointer(FuelModel), 
            PyMem_Malloc(300 * cython.sizeof(FuelModel)))
        if not self.fuel_models_arr:
            raise MemoryError()
        for f in fm.fuel_model_table.values():
            fm_number: pyidx = f["number"]
            fuel_model: FuelModel = fm_struct(f)
            self.fuel_models_arr[fm_number] = fuel_model
        
    @cy.ccall
    def get_fm_struct(self, fm_number: pyidx) -> FuelModel:
        fuel_models_arr: cy.pointer(FuelModel) = self.fuel_models_arr
        fuel_model: FuelModel = fuel_models_arr[fm_number]
        return fuel_model
    

    def __dealloc__(self):
        PyMem_Free(self.fuel_models_arr)  # no-op if self.data is NULL

@cy.ccall
def make_SpreadInputs(cube_resolution, space_time_cubes: dict) -> SpreadInputs:
    stc: SpreadInputs = SpreadInputs(cube_resolution,
        space_time_cubes["slope"],
        space_time_cubes["aspect"],
        space_time_cubes["fuel_model"],
        space_time_cubes["canopy_cover"],
        space_time_cubes["canopy_height"],
        space_time_cubes["canopy_base_height"],
        space_time_cubes["canopy_bulk_density"],
        space_time_cubes["wind_speed_10m"],
        space_time_cubes["upwind_direction"],
        space_time_cubes["fuel_moisture_dead_1hr"],
        space_time_cubes["fuel_moisture_dead_10hr"],
        space_time_cubes["fuel_moisture_dead_100hr"],
        space_time_cubes["fuel_moisture_live_herbaceous"],
        space_time_cubes["fuel_moisture_live_woody"],
        space_time_cubes["foliar_moisture"],
        space_time_cubes.get("temperature"),
        space_time_cubes.get("fuel_spread_adjustment"),
        space_time_cubes.get("weather_spread_adjustment")
    )
    return stc


@cy.profile(False)
@cy.cfunc
def lookup_space_time_cube_float32(space_time_cube: ISpaceTimeCube, space_time_coordinate: coord_tyx) -> cy.float:
    t: pyidx = space_time_coordinate[0]
    y: pyidx = space_time_coordinate[1]
    x: pyidx = space_time_coordinate[2]
    return space_time_cube.get_float(t, y, x)


CellInputs = cy.struct(
    slope                         = cy.float,
    aspect                        = cy.float,
    fuel_model_number             = cy.float,
    canopy_cover                  = cy.float,
    canopy_height                 = cy.float,
    canopy_base_height            = cy.float,
    canopy_bulk_density           = cy.float,
    wind_speed_10m                = cy.float,
    upwind_direction              = cy.float,
    fuel_moisture_dead_1hr        = cy.float,
    fuel_moisture_dead_10hr       = cy.float,
    fuel_moisture_dead_100hr      = cy.float,
    fuel_moisture_live_herbaceous = cy.float,
    fuel_moisture_live_woody      = cy.float,
    foliar_moisture               = cy.float,
    fuel_spread_adjustment        = cy.float,
    weather_spread_adjustment     = cy.float
)


@cy.cfunc
@cy.profile(True)
def lookup_cell_inputs(space_time_cubes: SpreadInputs, tyx: coord_tyx) -> CellInputs:

    (t, y, x) = tyx
    inputs: SpreadInputs = space_time_cubes
    # Topography, Fuel Model, and Vegetation
    slope               : cy.float = lookup_space_time_cube_float32(inputs.slope, tyx)               # rise/run
    aspect              : cy.float = lookup_space_time_cube_float32(inputs.aspect, tyx)              # degrees clockwise from North
    fuel_model_number   : cy.float = lookup_space_time_cube_float32(inputs.fuel_model, tyx)          # integer index in fm.fuel_model_table
    canopy_cover        : cy.float = lookup_space_time_cube_float32(inputs.canopy_cover, tyx)        # 0-1
    canopy_height       : cy.float = lookup_space_time_cube_float32(inputs.canopy_height, tyx)       # m
    canopy_base_height  : cy.float = lookup_space_time_cube_float32(inputs.canopy_base_height, tyx)  # m
    canopy_bulk_density : cy.float = lookup_space_time_cube_float32(inputs.canopy_bulk_density, tyx) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                : cy.float = lookup_space_time_cube_float32(inputs.wind_speed_10m, tyx)                # km/hr
    upwind_direction              : cy.float = lookup_space_time_cube_float32(inputs.upwind_direction, tyx)              # degrees clockwise from North
    fuel_moisture_dead_1hr        : cy.float = lookup_space_time_cube_float32(inputs.fuel_moisture_dead_1hr, tyx)        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       : cy.float = lookup_space_time_cube_float32(inputs.fuel_moisture_dead_10hr, tyx)       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      : cy.float = lookup_space_time_cube_float32(inputs.fuel_moisture_dead_100hr, tyx)      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous : cy.float = lookup_space_time_cube_float32(inputs.fuel_moisture_live_herbaceous, tyx) # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      : cy.float = lookup_space_time_cube_float32(inputs.fuel_moisture_live_woody, tyx)      # kg moisture/kg ovendry weight
    foliar_moisture               : cy.float = lookup_space_time_cube_float32(inputs.foliar_moisture, tyx)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    : cy.float = (lookup_space_time_cube_float32(inputs.fuel_spread_adjustment, tyx) 
                                 #if "fuel_spread_adjustment" in space_time_cubes
                                 if inputs.fuel_spread_adjustment is not None
                                 else 1.0)                                         # float >= 0.0
    weather_spread_adjustment : cy.float = (lookup_space_time_cube_float32(inputs.weather_spread_adjustment, tyx) 
                                 #if "weather_spread_adjustment" in space_time_cubes
                                 if inputs.weather_spread_adjustment is not None
                                 else 1.0)                                         # float >= 0.0

    return CellInputs(
        slope                          = slope,
        aspect                         = aspect,
        fuel_model_number              = fuel_model_number,
        canopy_cover                   = canopy_cover,
        canopy_height                  = canopy_height,
        canopy_base_height             = canopy_base_height,
        canopy_bulk_density            = canopy_bulk_density,
        wind_speed_10m                 = wind_speed_10m,
        upwind_direction               = upwind_direction,
        fuel_moisture_dead_1hr         = fuel_moisture_dead_1hr,
        fuel_moisture_dead_10hr        = fuel_moisture_dead_10hr,
        fuel_moisture_dead_100hr       = fuel_moisture_dead_100hr,
        fuel_moisture_live_herbaceous  = fuel_moisture_live_herbaceous,
        fuel_moisture_live_woody       = fuel_moisture_live_woody,
        foliar_moisture                = foliar_moisture,
        fuel_spread_adjustment         = fuel_spread_adjustment,
        weather_spread_adjustment      = weather_spread_adjustment
    )

@cy.ccall
@cy.cdivision(True)
@cy.exceptval(check=False)
def unburned_SpreadBehavior(elevation_gradient: vec_xy, dphi_st: vec_xyz) -> SpreadBehavior:
    # NOTE we're only going through these annoying calculations because we are required to return a spread_direction unit vector, which is of questionable value.
    # IMPROVEMENT We wouldn't have to go through this trouble if we simply returned a Cartesian speed vector instead, which would play more nicely with the rest of the code.
    spread_direction: vec_xyz
    if dphi_st != (0, 0, 0): # IMPROVEMENT numerical stability
        spread_direction = scale_3d(1./ vector_magnitude_3d(dphi_st), dphi_st)
    elif elevation_gradient != (0, 0): # IMPROVEMENT numerical stability
        slope_vector_3d: vec_xyz = to_slope_plane(elevation_gradient, elevation_gradient)
        spread_direction = as_unit_vector_3d(slope_vector_3d)
    else:
        spread_direction = (0., 1., 0.) # default: North
    #============================================================================================
    # Return zero surface fire behavior
    #============================================================================================
    return SpreadBehavior(
        dphi_dt           = 0.0,
        fire_type         = fire_type_unburned,
        spread_rate       = 0.0,
        spread_direction  = spread_direction,
        fireline_intensity= 0.0,
        flame_length      = 0.0,
    )

@cy.profile(True)
@cy.ccall
@cy.cdivision(True)
def burn_cell_toward_phi_gradient(space_time_cubes: object, 
                                  space_time_coordinate: coord_tyx, 
                                  phi_gradient_xy: vec_xy, 
                                  use_wind_limit: cy.bint = True,
                                  surface_lw_ratio_model: object = "rothermel",
                                  crown_max_lw_ratio: cy.float = 1e10 # FIXME None not allowed as float
                                  ) -> SpreadBehavior:
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
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
    - space_time_coordinate        :: (t,y,x)
    - phi_gradient_xy              :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)

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
    tyx: coord_tyx = (t, y, x)

    #================================================================================================
    # Unpack the space_time_cubes dictionary
    #================================================================================================

    inputs: SpreadInputs = space_time_cubes
    cell_inputs: CellInputs = lookup_cell_inputs(inputs, tyx)
    # Topography, Fuel Model, and Vegetation
    slope               : cy.float = cell_inputs.slope               # rise/run
    aspect              : cy.float = cell_inputs.aspect              # degrees clockwise from North
    fuel_model_number   : cy.float = cell_inputs.fuel_model_number   # integer index in fm.fuel_model_table
    canopy_cover        : cy.float = cell_inputs.canopy_cover        # 0-1
    canopy_height       : cy.float = cell_inputs.canopy_height       # m
    canopy_base_height  : cy.float = cell_inputs.canopy_base_height  # m
    canopy_bulk_density : cy.float = cell_inputs.canopy_bulk_density # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                : cy.float = cell_inputs.wind_speed_10m                # km/hr
    upwind_direction              : cy.float = cell_inputs.upwind_direction              # degrees clockwise from North
    fuel_moisture_dead_1hr        : cy.float = cell_inputs.fuel_moisture_dead_1hr        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       : cy.float = cell_inputs.fuel_moisture_dead_10hr       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      : cy.float = cell_inputs.fuel_moisture_dead_100hr      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous : cy.float = cell_inputs.fuel_moisture_live_herbaceous # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      : cy.float = cell_inputs.fuel_moisture_live_woody      # kg moisture/kg ovendry weight
    foliar_moisture               : cy.float = cell_inputs.foliar_moisture               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    : cy.float = cell_inputs.fuel_spread_adjustment    # float >= 0.0
    weather_spread_adjustment : cy.float = cell_inputs.weather_spread_adjustment # float >= 0.0
    spread_rate_adjustment    : cy.float = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

    #================================================================================================
    # Calculate the elevation gradient
    #================================================================================================

    elevation_gradient: vec_xy = calc_elevation_gradient(slope, aspect)

    #============================================================================================
    # Project the horizontal phi gradient onto the slope-tangential plane
    #============================================================================================

    phi_gradient: vec_xyz = calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient)

    #================================================================================================
    # Calculate the magnitude of the phi gradient
    #================================================================================================

    phi_magnitude: cy.float = vector_magnitude_3d(phi_gradient) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burnable
    #================================================================================================

    fm_number: pyidx = cy.cast(pyidx, fuel_model_number)
    #fm_struct: FuelModel = fuel_model_structs[fuel_model_number]
    fm_struct: FuelModel = inputs.get_fm_struct(fm_number)
    

    if not (phi_magnitude > 0.0 and fm_struct.burnable):
        return unburned_SpreadBehavior(elevation_gradient, phi_gradient)
    else:
        # Cell is on the fire perimeter and contains a burnable fuel model

        #============================================================================================
        # Compute derived parameters
        #============================================================================================

        M_f: fclaarr = (
            fuel_moisture_dead_1hr,
            fuel_moisture_dead_10hr,
            fuel_moisture_dead_100hr,
            0.0, # fuel_moisture_dead_herbaceous
            fuel_moisture_live_herbaceous,
            fuel_moisture_live_woody)
        fuel_bed_depth               : cy.float = fm_struct.delta                     # ft
        heat_of_combustion           : cy.float = Btu_lb_to_kJ_kg(fm_struct.h[0])     # kJ/kg
        estimated_fine_fuel_moisture : cy.float = fuel_moisture_dead_1hr              # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft: cy.float = wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min: cy.float = km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        # FIXME this is getting called as a Python function
        midflame_wind_speed: cy.float = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,  # m/min
                                                                    fuel_bed_depth,         # ft
                                                                    m_to_ft(canopy_height), # ft
                                                                    canopy_cover)           # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        #moisturized_fuel_model0 = fm.moisturize(fuel_model, fuel_moisture) # FIXME clean
        mfm: FuelModel = moisturize(fm_struct, M_f)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min: sf.FireBehaviorMin = sf.calc_surface_fire_behavior_no_wind_no_slope(mfm, spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max: sf.FireBehaviorMax = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                             midflame_wind_speed,
                                                             upwind_direction,
                                                             slope,
                                                             aspect,
                                                             use_wind_limit,
                                                             surface_lw_ratio_model)

        #============================================================================================
        # Calculate surface fire behavior normal to the fire perimeter
        #============================================================================================

        sfn: SpreadBehavior = calc_fireline_normal_behavior(surface_fire_max, phi_gradient)

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(sfn.fireline_intensity,
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            crown_fire_max: sf.FireBehaviorMax = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                             canopy_bulk_density, heat_of_combustion,
                                                             estimated_fine_fuel_moisture,
                                                             wind_speed_10m, upwind_direction,
                                                             slope, aspect, crown_max_lw_ratio)

            #========================================================================================
            # Calculate crown fire behavior normal to the fire perimeter
            #========================================================================================

            cfn: SpreadBehavior = calc_fireline_normal_behavior(crown_fire_max, phi_gradient)

            #========================================================================================
            # Calculate combined fire behavior normal to the fire perimeter
            #========================================================================================

            combined_fire_normal = cf.calc_combined_fire_behavior(sfn, cfn)
            
            #========================================================================================
            # Return the combined fire behavior normal to the fire perimeter
            #========================================================================================

            return combined_fire_normal

        else:

            #========================================================================================
            # Return the surface fire behavior normal to the fire perimeter
            #========================================================================================

            return sfn
# burn-cell-toward-phi-gradient ends here
# [[file:../../org/pyretechnics.org::phi-field-perimeter-tracking][phi-field-perimeter-tracking]]
@cy.profile(False)
@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def opposite_phi_signs(phi_matrix: cy.float[:,:], y1: pyidx, x1: pyidx, y2: pyidx, x2: pyidx) -> cy.bint:
    """
    Return True if the phi values at cells (x1,y1) and (x2,y2) have opposite signs.
    """
    return phi_matrix[y1, x1] * phi_matrix[y2, x2] < 0.0


# TODO: Is it faster to build up a list or a set?
# TODO: Should we store each frontier_cells entry as a coord_xy?
@cy.profile(True)
@cy.ccall
def identify_all_frontier_cells(phi_matrix: cy.float[:,:], rows: pyidx, cols: pyidx) -> set:
    """
    TODO: Add docstring
    """
    frontier_cells: set = set()
    y             : pyidx
    x             : pyidx
    for y in range(rows):
        for x in range(cols):
            # Compare (north, south, east, west) neighboring cell pairs for opposite phi signs
            north_y: pyidx = min(y+1, rows-1)
            south_y: pyidx = max(y-1, 0)
            east_x : pyidx = min(x+1, cols-1)
            west_x : pyidx = max(x-1, 0)
            if (opposite_phi_signs(phi_matrix, y, x, north_y, x) or
                opposite_phi_signs(phi_matrix, y, x, south_y, x) or
                opposite_phi_signs(phi_matrix, y, x, y, east_x) or
                opposite_phi_signs(phi_matrix, y, x, y, west_x)):
                frontier_cells.add((y, x))
    return frontier_cells

@cy.ccall
def is_frontier_cell(rows: pyidx, cols: pyidx, phi_matrix: cy.float[:,:], y: pyidx, x: pyidx) -> cy.bint:
    north_y: pyidx = min(y+1, rows-1)
    south_y: pyidx = max(y-1, 0)
    east_x : pyidx = min(x+1, cols-1)
    west_x : pyidx = max(x-1, 0)
    return (
        opposite_phi_signs(phi_matrix, y, x, north_y, x) or
        opposite_phi_signs(phi_matrix, y, x, south_y, x) or
        opposite_phi_signs(phi_matrix, y, x, y, east_x) or
        opposite_phi_signs(phi_matrix, y, x, y, west_x)
        )


# TODO: Is it faster to build up a list or a set?
# TODO: Should we store each frontier_cells entry as a coord_xy?
@cy.profile(True)
@cy.ccall
def identify_tracked_frontier_cells(phi_matrix: cy.float[:,:], fire_behavior_list: list, spot_ignited: list, rows: pyidx, cols: pyidx) -> set:
    """
    Resolves the set of frontier cells, returning a set of (y, x) coordinates.
    
    `fire_behavior_list` and `spot_ignited` are lists that cover the set of potential frontier cells.
    """
    frontier_cells: set = set()
    tcb: TrackedCellBehavior
    y: pyidx
    x: pyidx
    cell_index: object
    for tcb in fire_behavior_list:
        # Compare (north, south, east, west) neighboring cell pairs for opposite phi signs
        y = tcb.y
        x = tcb.x
        if is_frontier_cell(rows, cols, phi_matrix, y, x):
            cell_index = (y, x)
            frontier_cells.add(cell_index)
    for cell_index in spot_ignited:
        y = cell_index[0]
        x = cell_index[0]
        if is_frontier_cell(rows, cols, phi_matrix, y, x):
            frontier_cells.add(cell_index)
    return frontier_cells


@cy.profile(True)
@cy.ccall
def identify_tracked_cells(frontier_cells: set, buffer_width: pyidx, rows: pyidx, cols: pyidx) -> object:
    """
    TODO: Add docstring
    """
    tracked_cells: object = nbt.new_NarrowBandTracker(cols, rows)
    cell         : tuple
    for cell in frontier_cells:
        nbt.incr_square_around(tracked_cells, cell[0], cell[1], buffer_width)
    return tracked_cells


@cy.profile(True)
@cy.ccall
def update_tracked_cells(tracked_cells: object, frontier_cells_old: set, frontier_cells_new: set,
                         buffer_width: pyidx) -> object:
    """
    TODO: Add docstring
    """
    # Determine which frontier cells have been added or dropped
    frontier_cells_added  : set = frontier_cells_new.difference(frontier_cells_old)
    frontier_cells_dropped: set = frontier_cells_old.difference(frontier_cells_new)
    cell                  : tuple
    # Increment reference counters for all cells within buffer_width of the added frontier cells
    for cell in frontier_cells_added:
        nbt.incr_square_around(tracked_cells, cell[0], cell[1], buffer_width)
    # Decrement reference counters for all cells within buffer_width of the dropped frontier cells
    for cell in frontier_cells_dropped:
        nbt.decr_square_around(tracked_cells, cell[0], cell[1], buffer_width)
    # Return updated tracked cells
    return tracked_cells
# phi-field-perimeter-tracking ends here
# [[file:../../org/pyretechnics.org::spread-phi-field][spread-phi-field]]
# TODO: Move to pyretechnics.conversion
fire_type_codes = {
    "unburned"      : 0,
    "surface"       : 1,
    "passive_crown" : 2,
    "active_crown"  : 3,
}



@cy.cclass
class TrackedCellBehavior:
    spread_behavior: SpreadBehavior
    y: pyidx
    x: pyidx


@cy.inline
@cy.exceptval(check=False)
@cy.cfunc
def new_TrackedCellBehavior(
        spread_behavior: SpreadBehavior,
        y: pyidx,
        x: pyidx
        ) -> TrackedCellBehavior:
    t: TrackedCellBehavior = TrackedCellBehavior()
    t.y = y
    t.x = x
    t.spread_behavior = spread_behavior
    return t

@cy.ccall
def ignite_from_spotting(spot_ignitions: sortc.SortedDict, output_matrices, stop_time: cy.float) -> list:
    """
    Resolves the cells to be ignited by spotting in the current time step,
    returning them as a list of (y, x) tuples,
    and mutates `output_matrices` accordingly.
    """
    ignited: list = []
    if len(spot_ignitions) > 0:
        phi_matrix               : cy.float[:,:] = output_matrices["phi"]
        time_of_arrival_matrix   : cy.float[:,:] = output_matrices["time_of_arrival"]
        ignition_time: cy.float
        # https://grantjenks.com/docs/sortedcontainers/sorteddict.html
        n: pyidx = spot_ignitions.bisect_left(stop_time) # number of ignition_time values smaller than stop_time.
        _i: pyidx
        for _i in range(n):
            ignition_time, ignited_cells = spot_ignitions.popitem(index=0) # remove and return smallest ignition_time.
            for cell_index in ignited_cells:
                y: pyidx = cell_index[0]
                x: pyidx = cell_index[1]
                if phi_matrix[y,x] > 0.0: # Not burned yet
                    phi_matrix[y,x]             = -1.0
                    time_of_arrival_matrix[y,x] = ignition_time # FIXME: REVIEW Should I use stop_time instead?
                    ignited.append(cell_index)
                    # FIXME: I need to calculate and store the fire_behavior values for these cells
    return ignited

@cy.ccall
def reset_phi_star(
        fire_behavior_list: list,
        phi_star_matrix: cy.float[:,:],
        phi_matrix: cy.float[:,:],
        spot_ignited: list
    ) -> cy.void:
    y: pyidx
    x: pyidx
    for t in fire_behavior_list:
        tcb: TrackedCellBehavior = t
        y = tcb.y
        x = tcb.x
        phi_star_matrix[y,x] = phi_matrix[y,x]
    for cell_index in spot_ignited:
        y = cell_index[0]
        x = cell_index[1]
        phi_star_matrix[y,x] = phi_matrix[y,x]

@cy.ccall
def spot_from_burned_cell(
        spot_config,
        band_duration: cy.float,
        spatial_resolution: vec_xy,
        rows: pyidx, 
        cols: pyidx,
        stc: SpreadInputs,
        random_generator: BufferedRandGen,
        y: pyidx,
        x: pyidx,
        fb: SpreadBehavior,
        toa: cy.float,
        spot_ignitions: object
    ) -> cy.void:
    cell_horizontal_area_m2 : cy.float = spatial_resolution[0] * spatial_resolution[1]
    t_cast                  : pyidx    = int(toa // band_duration)
    space_time_coordinate   : coord_tyx = (t_cast, y, x)
    slope                   : cy.float = lookup_space_time_cube_float32(stc.slope, space_time_coordinate)
    aspect                  : cy.float = lookup_space_time_cube_float32(stc.aspect, space_time_coordinate)
    elevation_gradient      : vec_xy   = calc_elevation_gradient(slope, aspect)
    firebrands_per_unit_heat: cy.float = spot_config["firebrands_per_unit_heat"]
    expected_firebrand_count: cy.float = spot.expected_firebrand_production(fb,
                                                                            elevation_gradient,
                                                                            cell_horizontal_area_m2,
                                                                            firebrands_per_unit_heat)
    num_firebrands: pyidx = random_generator.next_poisson(expected_firebrand_count)
    if num_firebrands > 0:
        # OPTIM we might want to hold to the SpreadInputs and look these up in there.
        wind_speed_10m: cy.float = lookup_space_time_cube_float32(stc.wind_speed_10m, space_time_coordinate)
        upwind_direction: cy.float = lookup_space_time_cube_float32(stc.upwind_direction, space_time_coordinate)
        # FIXME optim first sample the number of firebrands (usually zero), then call this.
        new_ignitions: tuple[float, set]|None = spot.spread_firebrands(stc.fuel_model,
                                                                        stc.temperature,
                                                                        stc.fuel_moisture_dead_1hr,
                                                                        (rows, cols),
                                                                        spatial_resolution,
                                                                        space_time_coordinate,
                                                                        upwind_direction,
                                                                        wind_speed_10m,
                                                                        fb.fireline_intensity,
                                                                        fb.flame_length,
                                                                        toa,
                                                                        random_generator,
                                                                        num_firebrands,
                                                                        spot_config)
        if new_ignitions:
            ignition_time                      = new_ignitions[0]
            ignited_cells                      = new_ignitions[1]
            concurrent_ignited_cells: set|None = spot_ignitions.get(ignition_time)
            if concurrent_ignited_cells:
                spot_ignitions[ignition_time] = set.union(ignited_cells, concurrent_ignited_cells)
            else:
                spot_ignitions[ignition_time] = ignited_cells


@cy.profile(True)
# TODO: @cy.ccall
# TODO: @cy.exceptval(NULL)
# TODO: @cy.wraparound(False)
# TODO: @cy.boundscheck(False)
# TODO: Eliminate optional arguments to function
# TODO: Figure out how to type cube_resolution as vec_xyz
# TODO: Eliminate calls to np.dot
# TODO: See if we can use a C struct for the fire_behavior and fire_behavior_star dictionaries (not fire_behavior_dict)
# TODO: Perhaps cell_index should be a Python tuple rather than a C tuple?
# TODO: Replace fire_type_codes with as some kind of C type?
# TODO: Convert SpaceTimeCube to a @cclass with methods that take pyidx values
# TODO: Convert spot.spread_firebrands to a @cfunc
# TODO: Replace set.union with list concatenation? (Should ignited_cells be a list for speed?)
# TODO: Maybe ignition_time should be a Python float?
# TODO: Convert update_tracked_cells to a @cfunc
# TODO: Speed up burn_cell_toward_phi_gradient
# TODO: cimport the vec_xy and vec_xyz types to prevent typecasting when calling vector_magnitude_2d and dot_2d
# TODO: Turn off divide-by-zero checks
# TODO: Change for loops to use tracked_cells.keys() and sorted(spot_ignitions.keys())
# TODO: cimport numpy
@cy.ccall
def spread_fire_one_timestep(stc: SpreadInputs, output_matrices: dict, frontier_cells: set, tracked_cells: object,
                             cube_resolution: tuple, start_time: cy.float, max_timestep: cy.float,
                             max_cells_per_timestep: cy.float, use_wind_limit: bool|None = True,
                             surface_lw_ratio_model: str|None = "rothermel", crown_max_lw_ratio: float|None = None,
                             buffer_width: int|None = 3, spot_ignitions: dict|None = {}, spot_config: dict|None = None,
                             random_generator: object = None) -> dict:
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

    # Primitive defaults
    if use_wind_limit is None:
        use_wind_limit = True
    if surface_lw_ratio_model is None:
        surface_lw_ratio_model = "rothermel"
    if crown_max_lw_ratio is None:
        crown_max_lw_ratio = 1e10
    use_wind_limit: cy.bint = use_wind_limit
    surface_lw_ratio_model: object = surface_lw_ratio_model
    crown_max_lw_ratio: cy.float = crown_max_lw_ratio

    # Unpack output_matrices
    phi_matrix               : cy.float[:,:] = output_matrices["phi"]
    phi_star_matrix          : cy.float[:,:] = output_matrices["phi_star"]
    fire_type_matrix         : cy.uchar[:,:] = output_matrices["fire_type"]
    spread_rate_matrix       : cy.float[:,:] = output_matrices["spread_rate"]
    spread_direction_matrix  : cy.float[:,:] = output_matrices["spread_direction"]
    fireline_intensity_matrix: cy.float[:,:] = output_matrices["fireline_intensity"]
    flame_length_matrix      : cy.float[:,:] = output_matrices["flame_length"]
    time_of_arrival_matrix   : cy.float[:,:] = output_matrices["time_of_arrival"]

    # Extract simulation dimensions
    rows: pyidx = phi_matrix.shape[0]
    cols: pyidx = phi_matrix.shape[1]

    # Extract simulation resolution
    band_duration: cy.float = cube_resolution[0]
    cell_height  : cy.float = cube_resolution[1]
    cell_width   : cy.float = cube_resolution[2]
    spatial_resolution: vec_xy = (cell_height, cell_width)

    # Initialize max spread rates in the x and y dimensions to 0.0
    max_spread_rate_x: cy.float = 0.0
    max_spread_rate_y: cy.float = 0.0

    # Create an empty dictionary to store intermediate fire behavior values per cell
    fire_behavior_list: list = [] # INTRO a list recording the SpreadBehavior and tracked cells after the first pass.

    # First Runge-Kutta iteration.
    # Compute fire behavior values at start_time and identify the max spread rates in the x and y dimensions
    cell_index           : coord_yx
    y                    : pyidx
    x                    : pyidx
    space_time_coordinate: coord_tyx
    dphi_dt              : cy.float
    t0                   : pyidx = int(start_time // band_duration)
    tracked_cells_itr: nbt.TrackedCellsIterator = nbt.tracked_cells_iterator(tracked_cells)
    while tracked_cells_itr.has_next():
        cell_index = tracked_cells_itr.next_cell()
        # Unpack cell_index
        y = cell_index[0]
        x = cell_index[1]

        # Calculate phi gradient on the horizontal plane
        phi_gradient_xy : vec_xy   = calc_phi_gradient_approx(phi_matrix,
                                                              cell_width,
                                                              cell_height,
                                                              x,
                                                              y)
        phi_magnitude_xy: cy.float = vector_magnitude_2d(phi_gradient_xy)

        # Calculate the fire behavior normal to the fire front on the slope-tangential plane
        space_time_coordinate = (t0, y, x)
        fb: SpreadBehavior = burn_cell_toward_phi_gradient(stc,#space_time_cubes,
                                                           space_time_coordinate,
                                                           phi_gradient_xy,
                                                           use_wind_limit,
                                                           surface_lw_ratio_model,
                                                           crown_max_lw_ratio
                                                           )

        # Check whether cell has a positive phi magnitude
        if phi_magnitude_xy > 0.0:
            # Keep a running tally of the max horizontal spread rates in the x and y dimensions
            dphi_dt           : cy.float = fb.dphi_dt
            dphi_dx           : cy.float = phi_gradient_xy[0]
            dphi_dy           : cy.float = phi_gradient_xy[1]
            phi_magnitude_xy_2: cy.float = phi_magnitude_xy * phi_magnitude_xy
            spread_rate_x     : cy.float = -dphi_dt * dphi_dx / phi_magnitude_xy_2
            spread_rate_y     : cy.float = -dphi_dt * dphi_dy / phi_magnitude_xy_2
            max_spread_rate_x : cy.float = max(max_spread_rate_x, abs(spread_rate_x))
            max_spread_rate_y : cy.float = max(max_spread_rate_y, abs(spread_rate_y))

            # Integrate the Superbee flux limited phi gradient to make dphi_dt numerically stable
            phi_gradient_xy_limited: vec_xy = calc_phi_gradient(phi_matrix,
                                                                dphi_dx,
                                                                dphi_dy,
                                                                cell_width,
                                                                cell_height,
                                                                x,
                                                                y)
            dphi_dt_correction: cy.float = dot_2d(phi_gradient_xy, phi_gradient_xy_limited) / phi_magnitude_xy_2
            fb.dphi_dt = dphi_dt * dphi_dt_correction

        # Store fire behavior values for later use
        fire_behavior_list.append(new_TrackedCellBehavior(fb, y, x))

    # Calculate timestep using the CFL condition
    dt: cy.float
    if max_spread_rate_x == 0.0:
        if max_spread_rate_y == 0.0:
            dt = max_timestep
        else:
            dt = min(max_timestep, max_cells_per_timestep * cell_height / max_spread_rate_y)
    else:
        if max_spread_rate_y == 0.0:
            dt = min(max_timestep, max_cells_per_timestep * cell_width / max_spread_rate_x)
        else:
            dt = min(max_timestep, max_cells_per_timestep * min(cell_width / max_spread_rate_x,
                                                                cell_height / max_spread_rate_y))

    # Calculate the stop_time using this timestep
    stop_time: cy.float = start_time + dt

    # Update the tracked cell values in phi_star_matrix
    # FIXME factor out to function
    for t in fire_behavior_list:
        tcb: TrackedCellBehavior = t
        fb: SpreadBehavior = tcb.spread_behavior
        dphi_dt = fb.dphi_dt
        if dphi_dt != 0.0:
            phi_star_matrix[tcb.y,tcb.x] += dphi_dt * dt

    # Compute fire behavior values at stop_time and update the output_matrices
    t1           : pyidx = int(stop_time // band_duration)

    # 2nd Runge-Kutta iteration
    cell_horizontal_area_m2: cy.float = cube_resolution[1] * cube_resolution[2]
    ignition_time: float
    ignited_cells: set
    for t in fire_behavior_list:
        tcb: TrackedCellBehavior = t
        # Unpack cell_index
        y: pyidx = tcb.y
        x: pyidx = tcb.x

        # Calculate phi gradient on the horizontal plane
        phi_gradient_xy_star : vec_xy   = calc_phi_gradient_approx(phi_star_matrix,
                                                                   cell_width,
                                                                   cell_height,
                                                                   x,
                                                                   y)
        phi_magnitude_xy_star: cy.float = vector_magnitude_2d(phi_gradient_xy_star)

        # Calculate the fire behavior normal to the fire front on the slope-tangential plane
        space_time_coordinate    = (t1, y, x)
        fire_behavior_star: SpreadBehavior = burn_cell_toward_phi_gradient(stc,
                                                                 space_time_coordinate,
                                                                 phi_gradient_xy_star,
                                                                 use_wind_limit,
                                                                 surface_lw_ratio_model,
                                                                 crown_max_lw_ratio
                                                                 )

        # Check whether cell has a positive phi magnitude
        if phi_magnitude_xy_star > 0.0:
            # Integrate the Superbee flux limited phi gradient to make dphi_dt numerically stable
            dphi_dt_star                : cy.float = fire_behavior_star.dphi_dt
            dphi_dx_star                : cy.float = phi_gradient_xy_star[0]
            dphi_dy_star                : cy.float = phi_gradient_xy_star[1]
            phi_gradient_xy_star_limited: vec_xy   = calc_phi_gradient(phi_star_matrix,
                                                                       dphi_dx_star,
                                                                       dphi_dy_star,
                                                                       cell_width,
                                                                       cell_height,
                                                                       x,
                                                                       y)
            dphi_dt_star_correction: cy.float = (dot_2d(phi_gradient_xy_star, phi_gradient_xy_star_limited)
                                                 / (phi_magnitude_xy_star * phi_magnitude_xy_star))
            fire_behavior_star.dphi_dt = dphi_dt_star * dphi_dt_star_correction

        # Calculate the new phi value at stop_time as phi_next
        fb: SpreadBehavior = tcb.spread_behavior # NOTE we will use the fire behavior at the first Runge-Kutta iteration as output.
        dphi_dt_estimate1: cy.float = fb.dphi_dt
        dphi_dt_estimate2: cy.float = fire_behavior_star.dphi_dt
        # FIXME * 0.5
        dphi_dt_average  : cy.float = (dphi_dt_estimate1 + dphi_dt_estimate2) / 2.0
        if dphi_dt_average != 0.0: # FIXME why are we doing this unlikely check?
            phi     : cy.float = phi_matrix[y,x]
            phi_next: cy.float = phi + dphi_dt_average * dt

            # Update the tracked cell values in phi_matrix
            phi_matrix[y,x] = phi_next

            # Record fire behavior values in the output_matrices for cells that are burned in this timestep
            # NOTE: This records the fire behavior values at start_time and not at the time of arrival.
            if phi > 0.0 and phi_next <= 0.0: # This cell burned during this iteration.
                toa: cy.float = start_time + dt * phi / (phi - phi_next)
                fire_type_matrix[y,x]          = fb.fire_type
                spread_rate_matrix[y,x]        = fb.spread_rate
                spread_direction_matrix[y,x]   = spread_direction_vector_to_angle(fb.spread_direction)
                fireline_intensity_matrix[y,x] = fb.fireline_intensity
                flame_length_matrix[y,x]       = fb.flame_length
                time_of_arrival_matrix[y,x]    = toa

                # Cast firebrands, update firebrand_count_matrix, and update spot_ignitions
                # FIXME factor out to function
                if spot_config:
                    spot_from_burned_cell(
                        spot_config,
                        band_duration,
                        spatial_resolution,
                        rows, 
                        cols,
                        stc,
                        random_generator,
                        y,
                        x,
                        fb,
                        toa,
                        spot_ignitions
                        )

    spot_ignited = ignite_from_spotting(spot_ignitions, output_matrices, stop_time)

    # Save the new phi_matrix values in phi_star_matrix
    reset_phi_star(fire_behavior_list, phi_star_matrix, phi_matrix, spot_ignited)

    # Update the sets of frontier cells and tracked cells based on the updated phi matrix
    frontier_cells_new: set  = identify_tracked_frontier_cells(phi_matrix, fire_behavior_list, spot_ignited, rows, cols)
    tracked_cells_new : object = update_tracked_cells(tracked_cells, frontier_cells, frontier_cells_new, buffer_width)


    # Return the updated world state
    return {
        "simulation_time" : stop_time,
        "output_matrices" : output_matrices,
        "frontier_cells"  : frontier_cells_new,
        "tracked_cells"   : tracked_cells_new,
        "spot_ignitions"  : spot_ignitions,
        "random_generator": random_generator,
    }

@cy.ccall
@cy.inline
@cy.exceptval(check=False)
@cy.cdivision(True)
def compute_st_dphi_2(
    slp_dz: vec_xy, # Elevation gradient (dimensionless)
    dphi: vec_xy # 2D horizontal phi gradient (phi/m)
    ) -> cy.float:
    """
    Computes st_dphi_2, the squared norm (_2) of slope-tangential (st_) phi gradient (dphi).
    """
    dz0, dz1 = slp_dz
    dp0, dp1 = dphi
    p_z: cy.float = dp0*dz0 + dp1*dz1
    st_dphi_2: cy.float = dp0*dp0 + dp1*dp1 - (p_z * p_z)/(1.0 + dz0*dz0 + dz1*dz1)
    return st_dphi_2


PartialedEllWavelet = cy.struct( # INTRO Pre-computed coefficients to apply elliptical wavelet math as fast as possible (once the phi gradient information is available). See `dphi_dt_from_partialed_wavelet()``.
    Vh_3d = vec_xyz, # Heading spread rate vector (m/min).
    # The 'ewc_' prefix stands for 'elliptical wavelet coefficient', making the name easy to track for in the code.
    # I deliberately didn't choose descriptive names for these - I find that short, collision-free synthetic names are more useful to programmers here.
    ewc_A = cy.float, # Dimensionless coefficient (<= 0).
    ewc_B = cy.float, # Dimensionless coefficient (<= 0).
    ewc_C = cy.float, # Dimensionless coefficient (>= 0).
)


@cy.ccall
@cy.exceptval(check=False)
@cy.cdivision(True)
def prepare_partialed_wavelet(
        Vh_3d: vec_xyz, # Heading spread rate vector.
        Vf: cy.float, # Flanking spread rate.
        Vb: cy.float, # Backing spread rate.
        # NOTE accepting the (Vf, Vb) pair is more robust than just eccentricity (or equivalently LoW),
        # because not all fire models use elliptical wavelets that grow around their focus:
        # for example, ELMFIRE defaults to something else.
        ) -> PartialedEllWavelet:
    """
    Prepares ('partials') the elliptical wavelet calculation based on the given elliptical dimensions.
    The returned data is meant to be passed to function `dphi_dt_from_partialed_wavelet()`.
    """
    Vh: cy.float = vector_magnitude_3d(Vh_3d)
    Vh_inv: cy.float = (1.0/Vh) if Vh > 0 else 0.0
    LoW: cy.float = 0.5*(Vh + Vb)*Vh_inv # Length/Width ratio
    return PartialedEllWavelet(
        # NOTE if Vh = 0, all of the following will be 0, as they should.
        Vh_3d = Vh_3d, 
        ewc_A = -0.5 * (Vh - Vb) * Vh_inv, # INTRO
        ewc_B = -(Vf * Vh_inv), # INTRO
        ewc_C = (LoW*LoW - 1) # INTRO
        )


@cy.ccall
@cy.exceptval(check=False)
def zero_partialed_wavelet() -> PartialedEllWavelet:
    Vh_3d: vec_xyz = (0, 0, 0)
    return prepare_partialed_wavelet(Vh_3d, 0, 0)


@cy.ccall
@cy.cdivision(True)
@cy.exceptval(check=False)
def pw_from_FireBehaviorMax(fb_max: sf.FireBehaviorMax) -> PartialedEllWavelet:
    Vh: cy.float = fb_max.max_spread_rate
    Vh_3d: vec_xyz = scale_3d(Vh, fb_max.max_spread_direction)
    if Vh > 0:
        LoW: cy.float = fb_max.length_to_width_ratio
        ecc: cy.float = fb_max.eccentricity
        backing_adjustment: cy.float   = (1.0 - ecc) / (1.0 + ecc)
        Vb: cy.float  = Vh * backing_adjustment
        Vf: cy.float = (Vh + Vb) / (2.0 * LoW) # m/min
        return prepare_partialed_wavelet(Vh_3d, Vf, Vb)
    else:
        return zero_partialed_wavelet()



@cy.ccall
@cy.exceptval(check=False)
def dphi_dt_from_partialed_wavelet(
        pw: PartialedEllWavelet,
        dphi: vec_xy, # 2D horizontal gradient of phi.
        st_dphi_2: cy.float # INTRO Squared norm (_2) of slope-tangential (st_) phi gradient (dphi)
        ) -> cy.float:
    """
    Computes the dphi/dt (phi/min, <= 0) of one elliptical wavelet based on the spatial gradient of phi.
    """
    Vx, Vy, Vz = pw.Vh_3d
    Vh2 = (Vx*Vx + Vy*Vy + Vz*Vz) # Note: this could have been pre-computed too, but it's not clear that it would make things faster.
    dphi_dx, dphi_dy = dphi
    Delta: cy.float = (Vx*dphi_dx + Vy*dphi_dy) # Nontrivially, this is the dot-product (hence 'Delta') between the 3D heading spread rate vector and the 3D slope-tangential gradient of phi.
    dphi_dt: cy.float = (
        pw.ewc_A * Delta +
        pw.ewc_B * sqrt(
            Vh2 * st_dphi_2 +
            pw.ewc_C * (Delta*Delta)
            )
        )
    return dphi_dt

EllipticalInfo = cy.struct( # Pre-computed information required to compute dphi/dt, once the phi gradient is known. Derived from the surface and crowning wavelets.
    # NOTE the reason to make this a small struct stored in an array is efficiency - we want the CPU to have a low cache miss rate.
    cell_index = coord_yx,
    slp_dz = vec_xy, # Elevation gradient
    surfc_wavelet = PartialedEllWavelet,
    crowning_spread_rate = cy.float, # INTRO pre-computed critical threshold in surface spread rate at which crowning occurs.
    crown_wavelet = PartialedEllWavelet,
)
# NOTE a significant benefit of this architecture is that it's Rothermel-agnostic:
# EllipticalInfo could conceivably be implemented using variants of the Rothermel model.
# This can be valuable to give flexibility to users.

@cy.ccall
@cy.exceptval(check=False)
def dphi_dt_from_elliptical(ell_i: EllipticalInfo, dphi: vec_xy) -> cy.float: # NOTE this code has to be very fast!
    """
    Calculates the dphi/dt (a negative number in phi/min) of the combined surface and crown elliptical wavelets.

    The reason for computing and returning _only_ dphi/dt is efficiency:
    nothing else is needed in the front-propagating tight loop
    that iterates over tracked cells.
    """
    st_dphi_2: cy.float = compute_st_dphi_2(ell_i.slp_dz, dphi)
    surfc_dphi_dt: cy.float = dphi_dt_from_partialed_wavelet(ell_i.surfc_wavelet, dphi, st_dphi_2)
    csr: cy.float = ell_i.crowning_spread_rate
    csr_2: cy.float = csr * csr
    does_crown: cy.bint = (surfc_dphi_dt * surfc_dphi_dt) > csr_2 * st_dphi_2 # Logically equivalent to: (surface_spread_rate > crowning_spread_rate), but faster to compute and robust to zero phi gradient.
    if does_crown:
        crown_dphi_dt: cy.float = dphi_dt_from_partialed_wavelet(ell_i.crown_wavelet, dphi, st_dphi_2)
        return min(surfc_dphi_dt, crown_dphi_dt) # Note that dphi_dt <= 0
    else:
        return surfc_dphi_dt

p_CellInputs = cy.declare(pyidx, 17) # FIXME move

from cython.cimports.cython.view import array as cvarray # FIXME move

Pass1CellOutput = cy.struct( # INTRO some data saved during the 1st Runge-Kutta pass.
    cell_index = coord_yx,
    dphi = vec_xy,
    dphi_dt_0 = cy.float, # FIXME rename and be explicit about flux-limiting
    )

@cy.cclass
class TrackedCellsArrays:
    """
    Arrays used as on-heap supporting data structures during spread.
    """
    _arrays_length: pyidx # power of 2, double each time there is a dynamic resizing.
    n_tracked_cells: pyidx


    # FIXME timestamps that say when the data was last updated for each data column.

    # These are arrays
    float_inputs: cy.float[:, :] # Shape: (n_tracked_cells, p)
    ell_info: cy.pointer(EllipticalInfo) # Array of structs (needs to be iterated over very efficiently).
    pass1outputs: cy.pointer(Pass1CellOutput) # Per-cell data produced by the 1st Runge-Kutta pass.

    def __cinit__(self, _arrays_length: pyidx = 256):
        self._arrays_length = _arrays_length
        self.n_tracked_cells = 0
        self.float_inputs = cvarray(shape=(self._arrays_length, p_CellInputs), itemsize=cy.sizeof(cy.float), allocate_buffer=True, format='i')
        self.ell_info = cy.cast(cy.pointer(EllipticalInfo), PyMem_Malloc(self._arrays_length * cy.sizeof(EllipticalInfo)))
        self.pass1outputs = cy.cast(cy.pointer(Pass1CellOutput), PyMem_Malloc(self._arrays_length * cy.sizeof(Pass1CellOutput)))
        if not self.ell_info:
            raise MemoryError()


    @cy.ccall
    def reset_size(self, n_tracked_cells: pyidx) -> cy.void:
        """
        Ensures that this can hold at least `n_tracked_cells`, resizing the internal arrays if necessary.
        Also updates `self.n_tracked_cells`.
        This can erase any data present in this object: callers must make sure this information is no longer needed.
        """
        self.n_tracked_cells = n_tracked_cells
        while self.n_tracked_cells > self._arrays_length:
            self._arrays_length *= 2
        PyMem_Free(self.ell_info)
        self.ell_info = cy.cast(cy.pointer(EllipticalInfo), PyMem_Malloc(self._arrays_length * cy.sizeof(EllipticalInfo)))
        if not self.ell_info:
            raise MemoryError()
        self.float_inputs = cvarray(shape=(self._arrays_length, p_CellInputs), itemsize=cy.sizeof(cy.float), allocate_buffer=True, format='i')
        self.pass1outputs = cy.cast(cy.pointer(Pass1CellOutput), PyMem_Malloc(self._arrays_length * cy.sizeof(Pass1CellOutput)))
        if not self.pass1outputs:
            raise MemoryError()

    def __dealloc__(self):
        PyMem_Free(self.ell_info)


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def compare_cell_indexes(c0: coord_yx, c1: coord_yx) -> cy.int:
    """
    Lexicographic comparison of (y, x) coordinates. Returns -1, 0 or 1.
    """
    if c0[0] < c1[0]:
        return -1
    if c0[0] > c1[0]:
        return 1
    # Now we know both have the same y
    if c0[1] < c1[1]:
        return -1
    if c0[1] > c1[1]:
        return 1
    return 0


@cy.ccall
@cy.exceptval(check=False)
@cy.boundscheck(False)
@cy.wraparound(False)
def copy_tracked_cell_data(i_old: pyidx, tca_old: TrackedCellsArrays, i_new: pyidx, tca_new: TrackedCellsArrays) -> cy.void:
    for k in range(p_CellInputs): # OPTIM maybe we should unroll that loop, or use a native array copy function?
        tca_new.float_inputs[i_new, k] = tca_old.float_inputs[i_old, k]
    tca_new.ell_info[i_new] = tca_new.ell_info[i_old]
    # NOTE tca_old.pass1outputs does not need to be copied over given how it will get used.


@cy.cclass
class FireBehaviorSettings:
    """
    A fast-access data structure for fire behavior parameters,
    avoiding to pass lots of arguments around.
    """
    buffer_width: pyidx
    max_cells_per_timestep: cy.float # CFL condition

    use_wind_limit: cy.bint
    surface_lw_ratio_model: object
    crown_max_lw_ratio: cy.float
    
    spot_config: object


    def __init__(self, 
                max_cells_per_timestep: cy.float, 
                use_wind_limit: bool|None = True,
                surface_lw_ratio_model: str|None = "rothermel", 
                crown_max_lw_ratio: float|None = None,
                buffer_width: int|None = 3, 
                spot_config: dict|None = None) -> dict:
        self.buffer_width = (buffer_width if buffer_width is not None else 3)
        self.max_cells_per_timestep = max_cells_per_timestep
        
        self.use_wind_limit = (use_wind_limit if use_wind_limit is not None else True)
        self.surface_lw_ratio_model = (surface_lw_ratio_model if surface_lw_ratio_model is not None else "rothermel")
        self.crown_max_lw_ratio = (crown_max_lw_ratio if crown_max_lw_ratio is not None else 1e10)
        
        self.spot_config = spot_config


@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def load_float_inputs_for_cell(
        stc: SpreadInputs,
        tyx: coord_tyx,
        # Where to write the data to:
        tca: TrackedCellsArrays,
        i: pyidx) -> cy.void: # FIXME maybe return the CellInputs struct instead?
    """
    Reads variables from input SpaceTimeCubes and saves them by mutating `tca.float_inputs`.
    """
    float_inputs: cy.float[:,:] = tca.float_inputs
    float_inputs[i,0] = lookup_space_time_cube_float32(stc.slope, tyx)               # rise/run
    float_inputs[i,1] = lookup_space_time_cube_float32(stc.aspect, tyx)              # degrees clockwise from North
    
    float_inputs[i,2] = lookup_space_time_cube_float32(stc.fuel_model, tyx)          # integer index in fm.fuel_model_table
    
    float_inputs[i,3] = lookup_space_time_cube_float32(stc.canopy_cover, tyx)        # 0-1
    float_inputs[i,4] = lookup_space_time_cube_float32(stc.canopy_height, tyx)       # m
    float_inputs[i,5] = lookup_space_time_cube_float32(stc.canopy_base_height, tyx)  # m
    float_inputs[i,6] = lookup_space_time_cube_float32(stc.canopy_bulk_density, tyx) # kg/m^3
    
    float_inputs[i,7] = lookup_space_time_cube_float32(stc.wind_speed_10m, tyx)                # km/hr
    float_inputs[i,8] = lookup_space_time_cube_float32(stc.upwind_direction, tyx)              # degrees clockwise from North
    
    float_inputs[i,9] = lookup_space_time_cube_float32(stc.fuel_moisture_dead_1hr, tyx)        # kg moisture/kg ovendry weight
    float_inputs[i,10] = lookup_space_time_cube_float32(stc.fuel_moisture_dead_10hr, tyx)       # kg moisture/kg ovendry weight
    float_inputs[i,11] = lookup_space_time_cube_float32(stc.fuel_moisture_dead_100hr, tyx)      # kg moisture/kg ovendry weight
    float_inputs[i,12] = lookup_space_time_cube_float32(stc.fuel_moisture_live_herbaceous, tyx) # kg moisture/kg ovendry weight
    float_inputs[i,13] = lookup_space_time_cube_float32(stc.fuel_moisture_live_woody, tyx)      # kg moisture/kg ovendry weight
    float_inputs[i,14] = lookup_space_time_cube_float32(stc.foliar_moisture, tyx)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    float_inputs[i,15] = (lookup_space_time_cube_float32(stc.fuel_spread_adjustment, tyx) 
                                 if stc.fuel_spread_adjustment is not None
                                 else 1.0)                                         # float >= 0.0
    float_inputs[i,16] = (lookup_space_time_cube_float32(stc.weather_spread_adjustment, tyx) 
                                 if stc.weather_spread_adjustment is not None
                                 else 1.0)       


@cy.ccall
@cy.exceptval(check=False)
@cy.wraparound(False)
@cy.boundscheck(False)
def load_saved_CellInputs(float_inputs: cy.float[:,:], i: pyidx) -> CellInputs:
    """
    Loads the CellInputs struct by reading the data saved in the the float_inputs array.
    """
    slope               : cy.float = float_inputs[i, 0]               # rise/run
    aspect              : cy.float = float_inputs[i, 1]              # degrees clockwise from North
    
    fuel_model_number   : cy.float = float_inputs[i, 2]   # integer index in fm.fuel_model_table
    
    canopy_cover        : cy.float = float_inputs[i, 3]        # 0-1
    canopy_height       : cy.float = float_inputs[i, 4]       # m
    canopy_base_height  : cy.float = float_inputs[i, 5]  # m
    canopy_bulk_density : cy.float = float_inputs[i, 6] # kg/m^3
    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                : cy.float = float_inputs[i, 7]                # km/hr
    upwind_direction              : cy.float = float_inputs[i, 8]              # degrees clockwise from North
    
    fuel_moisture_dead_1hr        : cy.float = float_inputs[i, 9]        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       : cy.float = float_inputs[i, 10]       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      : cy.float = float_inputs[i, 11]      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous : cy.float = float_inputs[i, 12] # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      : cy.float = float_inputs[i, 13]      # kg moisture/kg ovendry weight
    foliar_moisture               : cy.float = float_inputs[i, 14]               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    : cy.float = float_inputs[i, 15]    # float >= 0.0
    weather_spread_adjustment : cy.float = float_inputs[i, 16] # float >= 0.0

    cell_inputs: CellInputs = CellInputs(
        slope = slope,
        aspect = aspect,
        fuel_model_number = fuel_model_number,
        canopy_cover = canopy_cover,
        canopy_height = canopy_height,
        canopy_base_height = canopy_base_height,
        canopy_bulk_density = canopy_bulk_density,
        wind_speed_10m = wind_speed_10m,
        upwind_direction = upwind_direction,
        fuel_moisture_dead_1hr = fuel_moisture_dead_1hr,
        fuel_moisture_dead_10hr = fuel_moisture_dead_10hr,
        fuel_moisture_dead_100hr = fuel_moisture_dead_100hr,
        fuel_moisture_live_herbaceous = fuel_moisture_live_herbaceous,
        fuel_moisture_live_woody = fuel_moisture_live_woody,
        foliar_moisture = foliar_moisture,
        fuel_spread_adjustment = fuel_spread_adjustment,
        weather_spread_adjustment = weather_spread_adjustment,
    )
    return cell_inputs


@cy.ccall
def resolve_surface_max_behavior(
        fb_opts: FireBehaviorSettings,
        ci: CellInputs,
        fm_struct: FuelModel,
        elevation_gradient: vec_xy,
        ) -> sf.FireBehaviorMax:
    spread_rate_adjustment    : cy.float = ci.fuel_spread_adjustment * ci.weather_spread_adjustment # float >= 0.0
    
    #============================================================================================
    # No wind/no-slope behavior
    #============================================================================================

    M_f: fclaarr = (
        ci.fuel_moisture_dead_1hr,
        ci.fuel_moisture_dead_10hr,
        ci.fuel_moisture_dead_100hr,
        0.0, # fuel_moisture_dead_herbaceous
        ci.fuel_moisture_live_herbaceous,
        ci.fuel_moisture_live_woody)
    # FIXME factor this out into a function for reusability.
    fuel_bed_depth               : cy.float = fm_struct.delta                     # ft
    # Apply fuel moisture to fuel model
    mfm: FuelModel = moisturize(fm_struct, M_f)
    surface_fire_min: sf.FireBehaviorMin = sf.calc_surface_fire_behavior_no_wind_no_slope(mfm, spread_rate_adjustment)

    #============================================================================================
    # Max surface fire behavior
    #============================================================================================

    # Convert from 10m wind speed to 20ft wind speed
    wind_speed_20ft: cy.float = wind_speed_10m_to_wind_speed_20ft(ci.wind_speed_10m) # km/hr

    # Convert 20ft wind speed from km/hr to m/min
    wind_speed_20ft_m_min: cy.float = km_hr_to_m_min(wind_speed_20ft) # m/min

    # Convert from 20ft wind speed to midflame wind speed in m/min
    # FIXME this is getting called as a Python function
    midflame_wind_speed: cy.float = sf.calc_midflame_wind_speed(
                                        wind_speed_20ft_m_min,  # m/min
                                        fuel_bed_depth,         # ft
                                        m_to_ft(ci.canopy_height), # ft
                                        ci.canopy_cover)           # 0-1
    
    # Calculate surface fire behavior in the direction of maximum spread
    surface_fire_max: sf.FireBehaviorMax = sf.calc_surface_fire_behavior_max(
                                                surface_fire_min,
                                                midflame_wind_speed,
                                                ci.upwind_direction,
                                                ci.slope,
                                                ci.aspect,
                                                fb_opts.use_wind_limit,
                                                fb_opts.surface_lw_ratio_model)
    return surface_fire_max

@cy.ccall
def resolve_crown_max_behavior(
        fb_opts: FireBehaviorSettings,
        ci: CellInputs,
        fm_struct: FuelModel,
        elevation_gradient: vec_xy,
        ) -> sf.FireBehaviorMax:
    heat_of_combustion: cy.float = Btu_lb_to_kJ_kg(fm_struct.h[0]) # kJ/kg
    return cf.calc_crown_fire_behavior_max(
        ci.canopy_height, 
        ci.canopy_base_height,
        ci.canopy_bulk_density, 
        heat_of_combustion,
        ci.wind_speed_10m, 
        ci.upwind_direction,
        ci.slope, 
        ci.aspect, 
        fb_opts.crown_max_lw_ratio)


@cy.ccall
@cy.cdivision(True)
def crown_critical_spread_rate( # FIXME maybe move to spot_fire module.
        ci: CellInputs,
        surface_fire_max: sf.FireBehaviorMax,
        ) -> cy.float:
    surfc_max_fli: cy.float = surface_fire_max.max_fireline_intensity
    crown_fli: cy.float = cf.van_wagner_critical_fireline_intensity(ci.canopy_base_height, ci.foliar_moisture)
    crown_spread_rate: cy.float = (surface_fire_max.max_spread_rate * crown_fli / surfc_max_fli # NOTE this uses the fact that fireline intensity is proportional to spread rate.
                                   if surfc_max_fli > 0
                                   else 0)
    return crown_spread_rate


@cy.ccall
@cy.cdivision(True)
@cy.wraparound(False)
@cy.boundscheck(False)
def resolve_cell_elliptical_info(
        fb_opts: FireBehaviorSettings,
        tyx: coord_tyx,
        stc: SpreadInputs, # NOTE only used to call stc.get_fm_struct(fm_number). We're not really reading the rasters here.
        float_inputs: cy.float[:,:],
        i: pyidx, 
        ) -> EllipticalInfo:
    
    cell_index: coord_yx = (tyx[1], tyx[2])
    ci: CellInputs = load_saved_CellInputs(float_inputs, i)
    elevation_gradient: vec_xy = calc_elevation_gradient(ci.slope, ci.aspect)

    # Load the fuel model
    fm_number: pyidx = cy.cast(pyidx, ci.fuel_model_number)
    #fm_struct: FuelModel = fuel_model_structs[fuel_model_number]
    fm_struct: FuelModel = stc.get_fm_struct(fm_number)

    surfc_pw: PartialedEllWavelet
    crown_pw: PartialedEllWavelet
    crown_spread_rate: cy.float
    if not fm_struct.burnable:
        surfc_pw = zero_partialed_wavelet()
        crown_spread_rate = 1234.5 # arbitrary positive value - this threshold will never be reached.
        crown_pw = zero_partialed_wavelet()
    else:
        surface_fire_max: sf.FireBehaviorMax = resolve_surface_max_behavior(fb_opts, ci, fm_struct, elevation_gradient)
        surfc_pw = pw_from_FireBehaviorMax(surface_fire_max)
        crown_spread_rate: cy.float = crown_critical_spread_rate(ci, surface_fire_max)
        crown_fire_max: sf.FireBehaviorMax = resolve_crown_max_behavior(fb_opts, ci, fm_struct, elevation_gradient)
        crown_pw = pw_from_FireBehaviorMax(crown_fire_max)
        
        
    #============================================================================================
    # Build the struct
    #============================================================================================

    ell_i: EllipticalInfo = EllipticalInfo(
        cell_index = cell_index,
        slp_dz = elevation_gradient,
        surfc_wavelet = surfc_pw,
        crowning_spread_rate = crown_spread_rate,
        crown_wavelet = crown_pw
    )
    return ell_i


@cy.ccall
def resolve_combined_spread_behavior(
        stc: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        tyx: coord_tyx,
        dphi: vec_xy,
        ) -> SpreadBehavior:
    """
    Similar to resolve_cell_elliptical_info, but does a more exhaustive computation and returns the SpreadBehavior struct.
    """
    ci: CellInputs = lookup_cell_inputs(stc, tyx)
    elevation_gradient: vec_xy = calc_elevation_gradient(ci.slope, ci.aspect)
    dphi_st: vec_xyz = calc_phi_gradient_on_slope(dphi, elevation_gradient)
    phi_magnitude: cy.float = vector_magnitude_3d(dphi_st)
    # Load the fuel model
    fm_number: pyidx = cy.cast(pyidx, ci.fuel_model_number)
    #fm_struct: FuelModel = fuel_model_structs[fuel_model_number]
    fm_struct: FuelModel = stc.get_fm_struct(fm_number)
    if phi_magnitude == 0.0 or not fm_struct.burnable:
        return unburned_SpreadBehavior(elevation_gradient, dphi_st)
    else:
        surface_fire_max: sf.FireBehaviorMax = resolve_surface_max_behavior(fb_opts, ci, fm_struct, elevation_gradient)
        crown_spread_rate: cy.float = crown_critical_spread_rate(ci, surface_fire_max)
        sfn: SpreadBehavior = calc_fireline_normal_behavior(surface_fire_max, dphi_st)
        doesnt_crown: cy.bint = (sfn.spread_rate <= crown_spread_rate)
        if doesnt_crown:
            return sfn
        else:
            crown_fire_max: sf.FireBehaviorMax = resolve_crown_max_behavior(fb_opts, ci, fm_struct, elevation_gradient)
            cfn: SpreadBehavior = calc_fireline_normal_behavior(crown_fire_max, dphi_st)
            combined_fire_normal: SpreadBehavior = cf.calc_combined_fire_behavior(sfn, cfn)
            return combined_fire_normal
        



@cy.ccall
def load_tracked_cell_data(
        # How to resolve inputs:
        stc: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        tyx: coord_tyx,
        # Where to write the data to:
        tca: TrackedCellsArrays,
        i: pyidx
        ) -> cy.void:
    load_float_inputs_for_cell(stc, tyx, tca, i)
    ell_i: EllipticalInfo = resolve_cell_elliptical_info(fb_opts, tyx, stc, tca.float_inputs, i)
    tca.ell_info[i] = ell_i



@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def sync_tracked_cells_arrays(
        stc: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        t: pyidx,
        tracked_cells: nbt.NarrowBandTracker,
        n_tracked_cells: pyidx,
        tca_old: TrackedCellsArrays,
        tca_new: TrackedCellsArrays
        ) -> cy.void:
    """
    Mutates `tca_new` to contain the same cells as `tracked_cells`, 
    """
    tca_new.reset_size(n_tracked_cells)
    cell_new: coord_yx
    i_old: pyidx = 0
    i_new: pyidx = 0
    cell_old: coord_yx = (0, 0)
    exhausted_old: cy.bint = i_old < tca_old.n_tracked_cells
    if not(exhausted_old):
        cell_old = tca_old.ell_info[i_old].cell_index
    # This loop uses the fact that both tca_old and new_cells_itr are sorted consistently with compare_cell_indexes().
    new_cells_itr: nbt.TrackedCellsIterator = nbt.tracked_cells_iterator(tracked_cells)
    while new_cells_itr.has_next(): # OPTIM if needed, make the iteration even faster by avoiding the use of the iterator, which will make the code a bit more ugly and coupled.
        cell_new = new_cells_itr.next_cell()
        while not(exhausted_old) and compare_cell_indexes(cell_old, cell_new) < 0: # cell_old is no longer tracked; just move forward.
            i_old += 1
            exhausted_old = i_old < tca_old.n_tracked_cells
            if not(exhausted_old):
                cell_old = tca_old.ell_info[i_old].cell_index
        if not(exhausted_old) and (compare_cell_indexes(cell_old, cell_new) == 0): # cell_new was already tracked: copy the data.
            copy_tracked_cell_data(i_old, tca_old, i_new, tca_new)
        else: # cell_new was not in tca_old
            tyx: coord_tyx = (t, cell_new[0], cell_new[1])
            load_tracked_cell_data(stc, fb_opts, tyx, tca_new, i_new)
        i_new += 1


@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def runge_kutta_pass1(
        fb_opts: FireBehaviorSettings, 
        spatial_resolution: vec_xy,
        max_timestep: cy.float,
        phi: cy.float[:, :],
        tca: TrackedCellsArrays
        ) -> cy.float:
    """
    1st Runge-Kutta loop over elliptical dimensions, which:
    1. Resolves dt from the CFL condition
    2. Saves a Pass1CellOutput struct for each cell.

    Returns the resolved `dt` and mutates `tca.pass1outputs`.
    Reads only `tca.ell_info`.
    """
    ell_info: cy.pointer[EllipticalInfo] = tca.ell_info # FIXME
    pass1outputs: cy.pointer[Pass1CellOutput] = tca.pass1outputs
    dt: cy.float = max_timestep
    (dx, dy) = spatial_resolution
    for i in range(tca.n_tracked_cells):
        ell_i: EllipticalInfo = ell_info[i]
        cell_index: coord_yx = ell_i.cell_index # FIXME
        y: pyidx = cell_index[0]
        x: pyidx = cell_index[1]
        dphi: vec_xy = calc_phi_gradient_approx(phi, dx, dy, x, y)
        dphi_dx: cy.float = dphi[0]
        dphi_dy: cy.float = dphi[1]
        dphi_norm2: cy.float = (dphi_dx * dphi_dx) + (dphi_dy * dphi_dy)
        dphi_dt_flim: cy.float
        if dphi_norm2 > 0: # Most common case.
            dphi_flim: vec_xy = calc_phi_gradient(phi, dphi_dx, dphi_dy, dx, dy, x, y) # Flux-limited 2D gradient.
            dphi_dt: cy.float = dphi_dt_from_elliptical(ell_i, dphi)
            dphi_dt_correction: cy.float = dot_2d(dphi, dphi_flim) / dphi_norm2
            dphi_dt_flim = (dphi_dt * dphi_dt_correction)
        else:
            dphi_dt_flim = 0.0
        pass1outputs[i] = Pass1CellOutput( # INTRO
            cell_index = cell_index,
            dphi = dphi,
            dphi_dt_0 = dphi_dt_flim
            )
        dt = dt # FIXME update dt, using some notion of min_ux
    return dt


@cy.ccall
@cy.wraparound(False)
@cy.boundscheck(False)
def update_phi_star(
        n_tracked_cells: pyidx, 
        tca: TrackedCellsArrays,
        dt: cy.float, 
        phi: cy.float[:,:], 
        phs: cy.float[:,:]
        ) -> cy.void:
    """
    Mutates the phs ('phi star') matrix, by using the dt and dphi/dt computed in the 1st Runge-Kutta pass.
    To be called between Runge-Kutta passes.
    """
    pass1outputs: cy.pointer[Pass1CellOutput] = tca.pass1outputs
    for i in range(n_tracked_cells):
        pass1out: Pass1CellOutput = pass1outputs[i]
        cell_index: coord_yx = pass1out.cell_index
        y, x = cell_index
        dphi_dt_0i: cy.float = pass1out.dphi_dt_0
        phs[y, x] = phi[y, x] + (dt * dphi_dt_0i)


@cy.cclass
class BurnedCellInfo: # Using an Extension Type instead of a struct because it's convenient to store in Python data structures like lists and dicts.
    """
    This data structure simply records information about a burned cell.
    """
    cell_index: coord_yx
    toa: cy.float # Time Of Arrival
    dphi: vec_xy # Gradient of Phi field, indicative of front direction.
    from_spotting: cy.bint # Whether spotting is what caused the cell to ignite.


@cy.ccall
def new_BurnedCellInfo( # Fast constructor function
        cell_index: coord_yx,
        toa: cy.float, # Time Of Arrival
        dphi: vec_xy, # Gradient of Phi field, indicative of front direction.
        from_spotting: cy.bint
        ) -> BurnedCellInfo:
    ret: BurnedCellInfo = BurnedCellInfo()
    ret.cell_index = cell_index
    ret.toa = toa
    ret.dphi = dphi
    ret.from_spotting = from_spotting
    return ret


@cy.ccall
@cy.cdivision(True)
@cy.wraparound(False)
@cy.boundscheck(False)
def runge_kutta_pass2(
        fb_opts: FireBehaviorSettings,
        spatial_resolution: vec_xy,
        start_time: cy.float,
        dt: cy.float,
        tca: TrackedCellsArrays,
        phi: cy.float[:, :],
        phs: cy.float[:, :]
        ) -> list[BurnedCellInfo]:
    """
    2nd Runge-Kutta loop, which:
    1. Updates the phi matrix
    2. Stores the old phi value
    3. Identifies cells that have just burned and returns them in a list.

    Reads from `tca` and `phs`, and mutates `phi`.
    """
    (dx, dy) = spatial_resolution
    n_tracked_cells: pyidx = tca.n_tracked_cells
    ell: cy.pointer[EllipticalInfo] = tca.ell_info
    pass1outputs: cy.pointer[Pass1CellOutput] = tca.pass1outputs
    spread_burned_cells: list[BurnedCellInfo] = []
    i: pyidx
    for i in range(n_tracked_cells):
        ell_i: EllipticalInfo = ell[i]
        cell_index: coord_yx = ell_i.cell_index
        y: pyidx = cell_index[0]
        x: pyidx = cell_index[1]
        dphi: vec_xy = calc_phi_gradient_approx(phs, dx, dy, x, y)
        dphi_dx: cy.float = dphi[0]
        dphi_dy: cy.float = dphi[1]
        dphi_norm2: cy.float = (dphi_dx * dphi_dx) + (dphi_dy * dphi_dy)
        dphi_dt_1i: cy.float
        if dphi_norm2 > 0: # Most common case.
            dphi_flim: vec_xy = calc_phi_gradient(phs, dphi_dx, dphi_dy, dx, dy, x, y) # Flux-limited 2D gradient.
            dphi_dt_correction: cy.float = dot_2d(dphi, dphi_flim) / dphi_norm2
            dphi_dt: cy.float = dphi_dt_from_elliptical(ell_i, dphi)
            dphi_dt_1i = (dphi_dt * dphi_dt_correction)
        else:
            dphi_dt_1i = 0.0
        dphi_dt_0i: cy.float = pass1outputs[i].dphi_dt_0
        phi_old: cy.float = phi[y, x] # NOTE could be inferred from phi_star.
        phi_new: cy.float = phi_old + 0.5 * dt * (dphi_dt_0i + dphi_dt_1i)
        phi[y, x] = phi_new
        i_just_burned: cy.bint = (phi_old * phi_new) < 0 # Phi can only ever decrease, therefore if these are of opposite signs, the cell has just burned.
        if i_just_burned:
            spread_burned_cells.append(new_BurnedCellInfo(
                cell_index, 
                toa = start_time + dt * phi_old / (phi_old - phi_new),
                dphi = pass1outputs[i].dphi, # Using the 1st-pass phi gradient. FIXME to be consistent with the toa, we might want to average this with the dphi from the 2nd pass.
                from_spotting=False))
    return spread_burned_cells


@cy.ccall
@cy.cdivision(True)
@cy.wraparound(False)
@cy.boundscheck(False)
def process_spread_burned_cells(
        stc: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        output_matrices: dict,
        spot_ignitions: object,
        random_generator: BufferedRandGen,
        spread_burned_cells: list[BurnedCellInfo]
        ) -> dict:
    fire_type_matrix: cy.uchar[:,:] = output_matrices["fire_type"]
    spread_rate_matrix: cy.float[:,:] = output_matrices["spread_rate"]
    spread_direction_matrix: cy.float[:,:] = output_matrices["spread_direction"]
    fireline_intensity_matrix: cy.float[:,:] = output_matrices["fireline_intensity"]
    flame_length_matrix: cy.float[:,:] = output_matrices["flame_length"]
    time_of_arrival_matrix: cy.float[:,:] = output_matrices["time_of_arrival"]

    burned_cell: BurnedCellInfo
    for burned_cell in spread_burned_cells:
        toa: cy.float = burned_cell.toa
        t: pyidx = int(toa // stc.band_duration)
        cell_index: coord_yx = burned_cell.cell_index
        y, x = cell_index
        tyx: coord_tyx = (t, y, x)
        # Re-compute the spread behavior. It's OK to re-compute it because a cell burning is a relatively rare event.
        fb: SpreadBehavior = resolve_combined_spread_behavior(stc, fb_opts, tyx, burned_cell.dphi)
        # Write to outputs
        fire_type_matrix[y,x]          = fb.fire_type
        spread_rate_matrix[y,x]        = fb.spread_rate
        spread_direction_matrix[y,x]   = spread_direction_vector_to_angle(fb.spread_direction)
        fireline_intensity_matrix[y,x] = fb.fireline_intensity
        flame_length_matrix[y,x]       = fb.flame_length
        time_of_arrival_matrix[y,x]    = toa

        # Cast firebrands, update firebrand_count_matrix, and update spot_ignitions
        if fb_opts.spot_config: # FIXME
            spot_from_burned_cell(
                fb_opts.spot_config,
                # FIXME simplify this function now that we have these new props:
                stc.band_duration,
                stc.spatial_resolution,
                stc.rows, 
                stc.cols,
                stc,
                random_generator,
                y,
                x,
                fb,
                toa,
                spot_ignitions)



# @cy.ccall
# @cy.cdivision(True)
# def spread_once_timestep(
#         sim_state: dict,
#         stc: SpreadInputs,
#         fb_opts: FireBehaviorSettings,
#         max_timestep: cy.float,
#     ) -> dict:
#     output_matrices: dict = sim_state["output_matrices"]
#     spot_ignitions: object = sim_state["spot_ignitions"]
#     random_generator: BufferedRandGen = sim_state["random_generator"]
#     phi: cy.float[:,:] = output_matrices["phi"]
#     phs: cy.float[:,:] = output_matrices["phi_star"]

#     band_duration: cy.float = stc.band_duration
#     (dx, dy) = stc.spatial_resolution
#     # FIXME resolve
#     n_tracked_cells: pyidx

#     i: pyidx # Index of tracked cells over 1D arrays
#     # Insert missing tracked cells.
#     tracked_cells: nbt.NarrowBandTracker = sim_state["tracked_cells"]
#     tca_old: TrackedCellsArrays = sim_state["_tracked_cells_arrays_old"]
#     tca: TrackedCellsArrays = sim_state["_tracked_cells_arrays"]
#     start_time: cy.float = sim_state["simulation_time"]
#     t0: pyidx = int(start_time // band_duration)
#     t_load: pyidx = t0
#     sync_tracked_cells_arrays(stc, fb_opts, t_load, tracked_cells, n_tracked_cells, tca_old, tca)
    
#     needs_recompute: cy.bint = False
#     # Refresh inputs if needed. FIXME "last_load_time" for each inputs.

#     # OPTIM: if fuel moistures haven't changed, don't recomute the no-wind/no-slope behavior.
#     # This makes sense if the wind field is refreshed more frequently than fuel moistures.
    
#     # Re-compute elliptical dimensions if needed. FIXME
    
#     dt = runge_kutta_pass1(fb_opts, stc.spatial_resolution, max_timestep, phi, tca)

#     # Now that dt is known, update phi_star_matrix.
#     update_phi_star(n_tracked_cells, tca, dt, phi, phs)
#     stop_time: cy.float = start_time + dt
    
#     spread_burned_cells: list[BurnedCellInfo] = runge_kutta_pass2(fb_opts, stc.spatial_resolution, start_time, dt, tca, phi, phs)
    
#     # Side-effects of the burned cells (outputs etc.).
#     process_spread_burned_cells(stc, fb_opts, output_matrices, spot_ignitions, random_generator, spread_burned_cells)
    

#     spot_ignited: list[BurnedCellInfo] = ignite_from_spotting(spot_ignitions, output_matrices, stop_time) # FIXME re-implement to produce BurnedCellInfo

#     # Save the new phi_matrix values in phi_star_matrix
#     reset_phi_star(spread_burned_cells, phs, phi, spot_ignited) # FIXME

#     # Update the sets of frontier cells and tracked cells based on the updated phi matrix
#     frontier_cells_new: set  = identify_tracked_frontier_cells(phi_matrix, fire_behavior_list, spot_ignited, rows, cols)
#     tracked_cells_new : object = update_tracked_cells(tracked_cells, frontier_cells, frontier_cells_new, buffer_width)


#     # Return the updated world state
#     sim_state_new = {
#         "simulation_time" : stop_time,
#         "output_matrices" : output_matrices,
#         "frontier_cells"  : frontier_cells_new,
#         "tracked_cells"   : tracked_cells_new,
#         # Purposefully swapping the tracked_cells_arrays
#         "_tracked_cells_arrays": tca_old,
#         "_tracked_cells_arrays_old": tca,
#         "spot_ignitions"  : spot_ignitions,
#         "random_generator": random_generator,
#     }
#     return sim_state_new



@cy.profile(True)
def spread_fire_with_phi_field(space_time_cubes, output_matrices, cube_resolution, start_time,
                               max_duration=None, use_wind_limit=True, surface_lw_ratio_model="rothermel",
                               crown_max_lw_ratio=None, max_cells_per_timestep=0.4, buffer_width=3,
                               spot_ignitions={}, spot_config=None):
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
      - slope                         :: rise/run
      - aspect                        :: degrees clockwise from North
      - fuel_model                    :: integer index in fm.fuel_model_table
      - canopy_cover                  :: 0-1
      - canopy_height                 :: m
      - canopy_base_height            :: m
      - canopy_bulk_density           :: kg/m^3
      - temperature                   :: degrees Celsius (Optional: needed for spotting)
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
    - output_matrices              :: dictionary of 2D Numpy arrays whose spatial dimensions match the space_time_cubes
      - phi                           :: 2D float array of values in [-1,1]
      - fire_type                     :: 2D byte array (0=unburned, 1=surface, 2=passive_crown, 3=active_crown)
      - spread_rate                   :: 2D float array (m/min)
      - spread_direction              :: 2D float array (degrees clockwise from North)
      - fireline_intensity            :: 2D float array (kW/m)
      - flame_length                  :: 2D float array (m)
      - time_of_arrival               :: 2D float array (min)
      - firebrand_count               :: 2D integer array (number of firebrands) (Optional)
    - cube_resolution              :: tuple with these fields
      - band_duration                 :: minutes
      - cell_height                   :: meters
      - cell_width                    :: meters
    - start_time                   :: minutes (from the start of the space_time_cube's temporal origin)
    - max_duration                 :: minutes (Optional)
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)
    - max_cells_per_timestep       :: max number of cells the fire front can travel in one timestep (Optional)
    - buffer_width                 :: Chebyshev distance from frontier cells to include in tracked cells (Optional)
    - spot_ignitions               :: dictionary of (ignition_time -> ignited_cells) (Optional: needed for spotting)
    - spot_config                  :: dictionary of spotting parameters (Optional: needed for spotting)
      - random_seed                   :: seed for a numpy.random.Generator object
      - firebrands_per_unit_heat      :: firebrands/kJ
      - downwind_distance_mean        :: meters
      - fireline_intensity_exponent   :: downwind_distance_mean multiplier [I^fireline_intensity_exponent]
      - wind_speed_exponent           :: downwind_distance_mean multiplier [U^wind_speed_exponent]
      - downwind_variance_mean_ratio  :: meters^2 / meter [downwind_variance_mean_ratio = Var(X) / E(X)]
      - crosswind_distance_stdev      :: meters
      - decay_distance                :: meters

    return a dictionary with these keys:
    - stop_time            :: minutes
    - stop_condition       :: "max duration reached" or "no burnable cells"
    - output_matrices      :: dictionary of 2D Numpy arrays whose spatial dimensions match the space_time_cubes
      - phi                   :: 2D float array of values in [-1,1]
      - fire_type             :: 2D byte array (0=unburned, 1=surface, 2=passive_crown, 3=active_crown)
      - spread_rate           :: 2D float array (m/min)
      - spread_direction      :: 2D float array (degrees clockwise from North)
      - fireline_intensity    :: 2D float array (kW/m)
      - flame_length          :: 2D float array (m)
      - time_of_arrival       :: 2D float array (min)
      - firebrand_count       :: 2D integer array (number of firebrands) (only included when provided as an input)
    - spot_ignitions       :: dictionary of (ignition_time -> ignited_cells) (only included when spotting is used)
    """
    # Extract simulation dimensions
    (bands, rows, cols) = space_time_cubes["slope"].shape
    band_duration       = cube_resolution[0]
    cube_duration       = bands * band_duration

    # Ensure that all space_time_cubes have the same spatial resolution
    for cube in space_time_cubes.values():
        if cube.shape != (bands, rows, cols):
            raise ValueError("The space_time_cubes must all share the same spatial resolution.")
    stc: SpreadInputs = make_SpreadInputs(cube_resolution, space_time_cubes)
    fb_opts: FireBehaviorSettings = FireBehaviorSettings( # FIXME pass to functions
        max_cells_per_timestep=max_cells_per_timestep,
        use_wind_limit = use_wind_limit,
        surface_lw_ratio_model = surface_lw_ratio_model,
        crown_max_lw_ratio = crown_max_lw_ratio,
        buffer_width = buffer_width,
        spot_config = spot_config,
        )

    # Ensure that space_time_cubes and output_matrices have the same spatial resolution
    for layer in output_matrices.values():
        if layer.shape != (rows, cols):
            raise ValueError("The space_time_cubes and output_matrices must share the same spatial resolution.")

    # Calculate the max stop time
    max_stop_time = start_time + max_duration if max_duration else cube_duration

    # Ensure that start_time does not exceed the cube_duration
    if start_time > cube_duration:
        raise ValueError("The start_time exceeds the temporal limit of the space_time_cubes.")

    # Ensure that the max_stop_time does not exceed the cube_duration
    if max_stop_time > cube_duration:
        raise ValueError("The start_time + max_duration exceeds the temporal limit of the space_time_cubes.")

    # Identify the sets of frontier cells and tracked cells based on the phi matrix
    phi_matrix     = output_matrices["phi"]
    frontier_cells = identify_all_frontier_cells(phi_matrix, rows, cols)
    tracked_cells  = identify_tracked_cells(frontier_cells, buffer_width, rows, cols)

    # Make a copy of the phi matrix to use for intermediate calculations in each timestep
    output_matrices["phi_star"] = np.copy(phi_matrix)

    # Create a numpy.random.Generator object to produce random samples if spot_config is provided
    random_generator = BufferedRandGen(np.random.default_rng(seed=spot_config["random_seed"])) if spot_config else None

    # Ensure that spot_ignitions is initialized as a SortedDict
    spot_igns = sortc.SortedDict(spot_ignitions)

    # Spread the fire until an exit condition is reached
    # FIXME: I don't think the "no burnable cells" condition can ever be met currently.
    simulation_time = start_time
    while(simulation_time < max_stop_time and (nbt.nonempty_tracked_cells(tracked_cells) or len(spot_ignitions) > 0)):
        # Compute max_timestep based on the remaining time in the temporal band and simulation
        remaining_time_in_band       = band_duration - simulation_time % band_duration
        remaining_time_in_simulation = max_stop_time - simulation_time
        max_timestep                 = min(remaining_time_in_band, remaining_time_in_simulation)

        # Spread fire one timestep
        results = spread_fire_one_timestep(stc, output_matrices, frontier_cells, tracked_cells,
                                           cube_resolution, simulation_time, max_timestep, max_cells_per_timestep,
                                           use_wind_limit, surface_lw_ratio_model, crown_max_lw_ratio,
                                           buffer_width, spot_igns, spot_config, random_generator)

        # Reset spread inputs
        simulation_time  = results["simulation_time"]
        output_matrices  = results["output_matrices"]
        frontier_cells   = results["frontier_cells"]
        tracked_cells    = results["tracked_cells"]
        spot_ignitions   = results["spot_ignitions"]
        random_generator = results["random_generator"]

    # Remove the temporary copy of the phi matrix from output_matrices
    output_matrices.pop("phi_star")

    # Return the final simulation results
    return {**{
        "stop_time"      : simulation_time,
        "stop_condition" : "max duration reached" if nbt.nonempty_tracked_cells(tracked_cells) else "no burnable cells",
        "output_matrices": output_matrices,
    }, **({ # FIXME restore | operator when done testing on older Python.
        "spot_ignitions"  : spot_igns,
        "random_generator": random_generator,
    } if spot_config else {})}
# spread-phi-field ends here
