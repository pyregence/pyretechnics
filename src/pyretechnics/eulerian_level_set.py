# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients-approx][phi-field-spatial-gradients-approx]]
import cython as cy
import numpy as np


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_dphi_dx_approx(phi: cy.float[:,:], dx: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                        cols: cy.Py_ssize_t) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given the cell width dx.
    """
    east_x: cy.Py_ssize_t = x + 1
    west_x: cy.Py_ssize_t = x - 1
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
@cy.exceptval(-1.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_dphi_dy_approx(phi: cy.float[:,:], dy: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                        rows: cy.Py_ssize_t) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given the cell height dy.
    """
    north_y: cy.Py_ssize_t = y + 1
    south_y: cy.Py_ssize_t = y - 1
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


# TODO: Pass rows and cols and create returned array without np.asarray?
# TODO: @cy.ccall
def calc_phi_gradient_approx(phi: cy.float[:,:], dx: cy.float, dy: cy.float,
                             x: cy.Py_ssize_t, y: cy.Py_ssize_t) -> cy.float[:,:]:
    """
    Calculate the spatial gradient of the phi raster at grid cell (x,y)
    given the cell width dx and the cell height dy.
    """
    rows   : cy.Py_ssize_t = phi.shape[0]
    cols   : cy.Py_ssize_t = phi.shape[1]
    dphi_dx: cy.float      = calc_dphi_dx_approx(phi, dx, x, y, cols)
    dphi_dy: cy.float      = calc_dphi_dy_approx(phi, dy, x, y, rows)
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
import cython as cy


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
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
import cython as cy
import numpy as np


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
def calc_dphi_dx(phi: cy.float[:,:], u_x: cy.float, dx: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                 cols: cy.Py_ssize_t) -> cy.float:
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
@cy.exceptval(-1.0)
def calc_dphi_dy(phi: cy.float[:,:], u_y: cy.float, dy: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                 rows: cy.Py_ssize_t) -> cy.float:
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


# TODO: Pass rows and cols and create returned array without np.asarray?
# TODO: @cy.ccall
def calc_phi_gradient(phi: cy.float[:,:], u_x: cy.float, u_y: cy.float, dx: cy.float, dy: cy.float,
                      x: cy.Py_ssize_t, y: cy.Py_ssize_t) -> cy.float[:,:]:
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
    rows   : cy.Py_ssize_t = phi.shape[0]
    cols   : cy.Py_ssize_t = phi.shape[1]
    dphi_dx: cy.float      = calc_dphi_dx(phi, u_x, dx, x, y, cols)
    dphi_dy: cy.float      = calc_dphi_dy(phi, u_y, dy, x, y, rows)
    return np.asarray((dphi_dx, dphi_dy))
# phi-field-spatial-gradients ends here
# [[file:../../org/pyretechnics.org::phi-east][phi-east]]
import cython as cy


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_east(phi: cy.float[:,:], u_x: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                  cols: cy.Py_ssize_t) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_x  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - cols :: integer number of columns in the phi matrix
    """
    very_east_x: cy.Py_ssize_t = min(x+2, cols-1)
    east_x     : cy.Py_ssize_t = min(x+1, cols-1)
    west_x     : cy.Py_ssize_t = max(x-1, 0)

    dphi_loc: cy.float = phi[y][east_x] - phi[y][x]
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
import cython as cy


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_west(phi: cy.float[:,:], u_x: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                  cols: cy.Py_ssize_t) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the -x (east->west)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_x  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - cols :: integer number of columns in the phi matrix
    """
    east_x     : cy.Py_ssize_t = min(x+1, cols-1)
    west_x     : cy.Py_ssize_t = max(x-1, 0)
    very_west_x: cy.Py_ssize_t = max(x-2, 0)

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
import cython as cy


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_north(phi: cy.float[:,:], u_y: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                   rows: cy.Py_ssize_t) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_y  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - rows :: integer number of rows in the phi matrix
    """
    very_north_y: cy.Py_ssize_t = min(y+2, rows-1)
    north_y     : cy.Py_ssize_t = min(y+1, rows-1)
    south_y     : cy.Py_ssize_t = max(y-1, 0)

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
import cython as cy


@cy.profile(False)
@cy.cfunc
@cy.exceptval(-1.0)
@cy.wraparound(False)
@cy.boundscheck(False)
def calc_phi_south(phi: cy.float[:,:], u_y: cy.float, x: cy.Py_ssize_t, y: cy.Py_ssize_t,
                   rows: cy.Py_ssize_t) -> cy.float:
    """
    Calculate the spatial gradient of the phi raster in the -y (north->south)
    direction at grid cell (x,y) given:
    - phi  :: 2D float array of values in [-1,1]
    - u_y  :: m/min
    - x    :: integer column index in phi
    - y    :: integer row index in phi
    - rows :: integer number of rows in the phi matrix
    """
    north_y     : cy.Py_ssize_t = min(y+1, rows-1)
    south_y     : cy.Py_ssize_t = max(y-1, 0)
    very_south_y: cy.Py_ssize_t = max(y-2, 0)

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
from math import sqrt
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.surface_fire as sf
import pyretechnics.vector_utils as vu


# TODO: Move this to pyretechnics.vector_utils and use throughout the literate program
def calc_elevation_gradient(slope, aspect):
    """
    Returns the elevation gradient (dz_dx: rise/run, dz_dy: rise/run) given:
    - slope  :: rise/run
    - aspect :: degrees clockwise from North
    """
    return conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))


def calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient):
    """
    Return the gradient of phi projected onto the slope-tangential plane as a 3D (x,y,z) vector (in phi/m) given:
    - phi_gradient_xy    :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - elevation_gradient :: (dz_dx: rise/run, dz_dy: rise/run)
    """
    (dphi_dx, dphi_dy) = phi_gradient_xy
    phi_gradient_xyz   = np.asarray((dphi_dx, dphi_dy, 0.0))
    if vu.vector_magnitude(elevation_gradient) == 0.0:
        return phi_gradient_xyz
    else:
        slope_normal_vector = vu.get_slope_normal_vector(elevation_gradient) # (x,y,z) unit vector
        return phi_gradient_xyz - np.dot(phi_gradient_xyz, slope_normal_vector) * slope_normal_vector


# FIXME: Do I switch to cruz_passive_crown_fire_spread_rate() if the normal_spread_rate < critical_spread_rate?
#        Did I do this correctly in calc_crown_fire_behavior_in_direction?
def calc_fireline_normal_behavior(fire_behavior_max, phi_gradient):
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

    phi_magnitude = vu.vector_magnitude(phi_gradient) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burning
    #================================================================================================

    if (phi_magnitude == 0.0 or fire_behavior_max["max_spread_rate"] == 0.0):
        # This location is not on the fire perimeter and/or is not burning

        #================================================================================================
        # Set the spread direction to the phi gradient direction, max spread direction, upslope, or North
        #================================================================================================

        spread_direction = (np.asarray(phi_gradient) / phi_magnitude if phi_magnitude > 0.0
                            else fire_behavior_max["max_spread_direction"])

        #============================================================================================
        # Return zero surface/crown fire behavior
        #============================================================================================

        return {
            "dphi_dt"           : 0.0,
            "fire_type"         : "unburned",
            "spread_rate"       : 0.0,
            "spread_direction"  : spread_direction,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }

    else:
        # This location is on the fire perimeter and is burning

        #============================================================================================
        # Unpack the fire_behavior_max dictionary
        #============================================================================================

        heading_fire_type          = fire_behavior_max.get("max_fire_type", "surface")
        heading_spread_rate        = fire_behavior_max["max_spread_rate"]               # m/min
        heading_spread_direction   = fire_behavior_max["max_spread_direction"]          # (x,y,z) unit vector
        heading_spread_vector      = heading_spread_rate * heading_spread_direction     # (x,y,z) m/min vector
        heading_fireline_intensity = fire_behavior_max["max_fireline_intensity"]        # kW/m
        length_to_width_ratio      = fire_behavior_max["length_to_width_ratio"]         # unitless
        eccentricity               = fire_behavior_max["eccentricity"]                  # unitless
        critical_spread_rate       = fire_behavior_max.get("critical_spread_rate", 0.0) # m/min

        #============================================================================================
        # Calculate the backing and flanking fire spread rates
        #============================================================================================

        backing_adjustment   = (1.0 - eccentricity) / (1.0 + eccentricity)                                 # unitless
        backing_spread_rate  = heading_spread_rate * backing_adjustment                                    # m/min
        flanking_spread_rate = (heading_spread_rate + backing_spread_rate) / (2.0 * length_to_width_ratio) # m/min

        #============================================================================================
        # Calculate dphi/dt
        #============================================================================================

        A       = (heading_spread_rate - backing_spread_rate) / (2 * heading_spread_rate) # unitless
        B       = np.dot(heading_spread_vector, phi_gradient)                             # phi/min
        C       = flanking_spread_rate / heading_spread_rate                              # unitless
        D       = (heading_spread_rate * phi_magnitude) ** 2.0                            # (phi/min)^2
        E       = (length_to_width_ratio ** 2.0 - 1.0) * (B ** 2.0)                       # (phi/min)^2
        dphi_dt = -(A * B + C * sqrt(D + E))                                              # phi/min

        #============================================================================================
        # Calculate fire behavior normal to the fire perimeter
        #============================================================================================

        normal_spread_rate        = -dphi_dt / phi_magnitude                        # m/min
        normal_direction          = np.asarray(phi_gradient) / phi_magnitude        # (x,y,z) unit vector
        normal_adjustment         = normal_spread_rate / heading_spread_rate        # unitless
        normal_fireline_intensity = heading_fireline_intensity * normal_adjustment  # kW/m
        normal_flame_length       = sf.calc_flame_length(normal_fireline_intensity) # m
        normal_fire_type          = ("surface" if heading_fire_type == "surface"
                                     else "active_crown" if normal_spread_rate > critical_spread_rate
                                     else "passive_crown")

        #========================================================================================
        # Return the surface/crown fire behavior normal to the fire perimeter
        #========================================================================================

        return {
            "dphi_dt"           : dphi_dt,                   # phi/min
            "fire_type"         : normal_fire_type,          # surface, passive_crown, or active_crown
            "spread_rate"       : normal_spread_rate,        # m/min
            "spread_direction"  : normal_direction,          # (x,y,z) unit vector
            "fireline_intensity": normal_fireline_intensity, # kW/m
            "flame_length"      : normal_flame_length,       # m
        }
# calc-fireline-normal-behavior ends here
# [[file:../../org/pyretechnics.org::burn-cell-toward-phi-gradient][burn-cell-toward-phi-gradient]]
import numpy as np
import pyretechnics.conversion as conv
import pyretechnics.crown_fire as cf
import pyretechnics.fuel_models as fm
import pyretechnics.surface_fire as sf
import pyretechnics.vector_utils as vu


# TODO: Create a version of this function that runs efficiently over a space_time_region
def burn_cell_toward_phi_gradient(space_time_cubes, space_time_coordinate, phi_gradient_xy, use_wind_limit=True,
                                  surface_lw_ratio_model="rothermel", crown_max_lw_ratio=None):
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

    #================================================================================================
    # Unpack the space_time_cubes dictionary
    #================================================================================================

    # Topography, Fuel Model, and Vegetation
    slope               = space_time_cubes["slope"].get(t,y,x)               # rise/run
    aspect              = space_time_cubes["aspect"].get(t,y,x)              # degrees clockwise from North
    fuel_model_number   = space_time_cubes["fuel_model"].get(t,y,x)          # integer index in fm.fuel_model_table
    canopy_cover        = space_time_cubes["canopy_cover"].get(t,y,x)        # 0-1
    canopy_height       = space_time_cubes["canopy_height"].get(t,y,x)       # m
    canopy_base_height  = space_time_cubes["canopy_base_height"].get(t,y,x)  # m
    canopy_bulk_density = space_time_cubes["canopy_bulk_density"].get(t,y,x) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m                = space_time_cubes["wind_speed_10m"].get(t,y,x)                # km/hr
    upwind_direction              = space_time_cubes["upwind_direction"].get(t,y,x)              # degrees clockwise from North
    fuel_moisture_dead_1hr        = space_time_cubes["fuel_moisture_dead_1hr"].get(t,y,x)        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr       = space_time_cubes["fuel_moisture_dead_10hr"].get(t,y,x)       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr      = space_time_cubes["fuel_moisture_dead_100hr"].get(t,y,x)      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous = space_time_cubes["fuel_moisture_live_herbaceous"].get(t,y,x) # kg moisture/kg ovendry weight
    fuel_moisture_live_woody      = space_time_cubes["fuel_moisture_live_woody"].get(t,y,x)      # kg moisture/kg ovendry weight
    foliar_moisture               = space_time_cubes["foliar_moisture"].get(t,y,x)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment    = (space_time_cubes["fuel_spread_adjustment"].get(t,y,x)
                                 if "fuel_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    weather_spread_adjustment = (space_time_cubes["weather_spread_adjustment"].get(t,y,x)
                                 if "weather_spread_adjustment" in space_time_cubes
                                 else 1.0)                                         # float >= 0.0
    spread_rate_adjustment    = fuel_spread_adjustment * weather_spread_adjustment # float >= 0.0

    #================================================================================================
    # Calculate the elevation gradient
    #================================================================================================

    elevation_gradient = calc_elevation_gradient(slope, aspect)

    #============================================================================================
    # Project the horizontal phi gradient onto the slope-tangential plane
    #============================================================================================

    phi_gradient = calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient)

    #================================================================================================
    # Calculate the magnitude of the phi gradient
    #================================================================================================

    phi_magnitude = vu.vector_magnitude(phi_gradient) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burnable
    #================================================================================================

    fuel_model = fm.fuel_model_table.get(fuel_model_number)

    if not (phi_magnitude > 0.0 and fuel_model and fuel_model["burnable"]):
        # Cell is not on the fire perimeter and/or contains an unknown or non-burnable fuel model

        #================================================================================================
        # Set the spread direction to the phi gradient direction, upslope, or North
        #================================================================================================

        if phi_magnitude > 0.0:
            spread_direction = np.asarray(phi_gradient) / phi_magnitude
        elif slope > 0.0:
            slope_vector_3d  = vu.to_slope_plane(elevation_gradient, elevation_gradient)
            spread_direction = vu.as_unit_vector(slope_vector_3d)
        else:
            spread_direction = np.asarray((0,1,0)) # default: North

        #============================================================================================
        # Return zero surface fire behavior
        #============================================================================================

        return {
            "dphi_dt"           : 0.0,
            "fire_type"         : "unburned",
            "spread_rate"       : 0.0,
            "spread_direction"  : spread_direction,
            "fireline_intensity": 0.0,
            "flame_length"      : 0.0,
        }

    else:
        # Cell is on the fire perimeter and contains a burnable fuel model

        #============================================================================================
        # Compute derived parameters
        #============================================================================================

        fuel_moisture                = [fuel_moisture_dead_1hr,
                                        fuel_moisture_dead_10hr,
                                        fuel_moisture_dead_100hr,
                                        0.0, # fuel_moisture_dead_herbaceous
                                        fuel_moisture_live_herbaceous,
                                        fuel_moisture_live_woody]               # kg moisture/kg ovendry weight
        fuel_bed_depth               = fuel_model["delta"]                      # ft
        heat_of_combustion           = conv.Btu_lb_to_kJ_kg(fuel_model["h"][0]) # kJ/kg
        estimated_fine_fuel_moisture = fuel_moisture_dead_1hr                   # kg moisture/kg ovendry weight

        #============================================================================================
        # Calculate midflame wind speed
        #============================================================================================

        # Convert from 10m wind speed to 20ft wind speed
        wind_speed_20ft = conv.wind_speed_10m_to_wind_speed_20ft(wind_speed_10m) # km/hr

        # Convert 20ft wind speed from km/hr to m/min
        wind_speed_20ft_m_min = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

        # Convert from 20ft wind speed to midflame wind speed in m/min
        midflame_wind_speed = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,       # m/min
                                                          fuel_bed_depth,              # ft
                                                          conv.m_to_ft(canopy_height), # ft
                                                          canopy_cover)                # 0-1

        #============================================================================================
        # Calculate surface fire behavior in the direction of maximum spread
        #============================================================================================

        # Apply fuel moisture to fuel model
        moisturized_fuel_model = fm.moisturize(fuel_model, fuel_moisture)

        # TODO: Memoize calc_surface_fire_behavior_no_wind_no_slope
        # Calculate no-wind-no-slope surface fire behavior
        surface_fire_min = sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model,
                                                                          spread_rate_adjustment)

        # Calculate surface fire behavior in the direction of maximum spread
        surface_fire_max = sf.calc_surface_fire_behavior_max(surface_fire_min,
                                                             midflame_wind_speed,
                                                             upwind_direction,
                                                             slope,
                                                             aspect,
                                                             use_wind_limit,
                                                             surface_lw_ratio_model)

        #============================================================================================
        # Calculate surface fire behavior normal to the fire perimeter
        #============================================================================================

        surface_fire_normal = calc_fireline_normal_behavior(surface_fire_max, phi_gradient)

        #============================================================================================
        # Determine whether the surface fire transitions to a crown fire
        #============================================================================================

        if cf.van_wagner_crown_fire_initiation(surface_fire_normal["fireline_intensity"],
                                               canopy_cover,
                                               canopy_base_height,
                                               foliar_moisture):

            #========================================================================================
            # Calculate crown fire behavior in the direction of maximum spread
            #========================================================================================

            crown_fire_max = cf.calc_crown_fire_behavior_max(canopy_height, canopy_base_height,
                                                             canopy_bulk_density, heat_of_combustion,
                                                             estimated_fine_fuel_moisture,
                                                             wind_speed_10m, upwind_direction,
                                                             slope, aspect, crown_max_lw_ratio)

            #========================================================================================
            # Calculate crown fire behavior normal to the fire perimeter
            #========================================================================================

            crown_fire_normal = calc_fireline_normal_behavior(crown_fire_max, phi_gradient)

            #========================================================================================
            # Calculate combined fire behavior normal to the fire perimeter
            #========================================================================================

            combined_fire_normal = cf.calc_combined_fire_behavior(surface_fire_normal, crown_fire_normal)
            surface_dphi_dt      = surface_fire_normal["dphi_dt"]
            crown_dphi_dt        = crown_fire_normal["dphi_dt"]
            combined_dphi_dt     = surface_dphi_dt if abs(surface_dphi_dt) > abs(crown_dphi_dt) else crown_dphi_dt
            combined_fire_normal["dphi_dt"] = combined_dphi_dt

            #========================================================================================
            # Return the combined fire behavior normal to the fire perimeter
            #========================================================================================

            return combined_fire_normal

        else:

            #========================================================================================
            # Return the surface fire behavior normal to the fire perimeter
            #========================================================================================

            return surface_fire_normal
# burn-cell-toward-phi-gradient ends here
# [[file:../../org/pyretechnics.org::phi-field-perimeter-tracking][phi-field-perimeter-tracking]]
def opposite_phi_signs(phi_matrix, y1, x1, y2, x2):
    """
    TODO: Add docstring
    """
    return phi_matrix[y1, x1] * phi_matrix[y2, x2] < 0.0


def identify_frontier_cells(phi_matrix, tracked_cells=None):
    """
    TODO: Add docstring
    """
    (rows, cols)   = phi_matrix.shape
    frontier_cells = set()
    if tracked_cells:
        for (y, x) in tracked_cells:
            # Compare (north, south, east, west) neighboring cell pairs for opposite phi signs
            north_y = min(y+1, rows-1)
            south_y = max(y-1, 0)
            east_x  = min(x+1, cols-1)
            west_x  = max(x-1, 0)
            if (opposite_phi_signs(phi_matrix, y, x, north_y, x) or
                opposite_phi_signs(phi_matrix, y, x, south_y, x) or
                opposite_phi_signs(phi_matrix, y, x, y, east_x) or
                opposite_phi_signs(phi_matrix, y, x, y, west_x)):
                frontier_cells.add((y, x))
    else:
        for y in range(rows):
            for x in range(cols):
                # Compare (north, south, east, west) neighboring cell pairs for opposite phi signs
                north_y = min(y+1, rows-1)
                south_y = max(y-1, 0)
                east_x  = min(x+1, cols-1)
                west_x  = max(x-1, 0)
                if (opposite_phi_signs(phi_matrix, y, x, north_y, x) or
                    opposite_phi_signs(phi_matrix, y, x, south_y, x) or
                    opposite_phi_signs(phi_matrix, y, x, y, east_x) or
                    opposite_phi_signs(phi_matrix, y, x, y, west_x)):
                    frontier_cells.add((y, x))
    return frontier_cells


def project_buffer(cell, buffer_width, rows, cols):
    """
    TODO: Add docstring
    """
    (y, x)         = cell
    buffer_range_y = range(max(0, y - buffer_width), min(rows, y + buffer_width + 1))
    buffer_range_x = range(max(0, x - buffer_width), min(cols, x + buffer_width + 1))
    return [(y_, x_)
            for y_ in buffer_range_y
            for x_ in buffer_range_x]


def identify_tracked_cells(frontier_cells, buffer_width, matrix_shape):
    """
    TODO: Add docstring
    """
    (rows, cols)  = matrix_shape
    tracked_cells = {}
    for cell in frontier_cells:
        for buffer_cell in project_buffer(cell, buffer_width, rows, cols):
            tracked_cells[buffer_cell] = tracked_cells.get(buffer_cell, 0) + 1
    return tracked_cells


def update_tracked_cells(tracked_cells, frontier_cells_old, frontier_cells_new, buffer_width, matrix_shape):
    """
    TODO: Add docstring
    """
    # Determine which frontier cells have been added or dropped
    (rows, cols)           = matrix_shape
    frontier_cells_added   = frontier_cells_new.difference(frontier_cells_old)
    frontier_cells_dropped = frontier_cells_old.difference(frontier_cells_new)
    # Increment reference counters for all cells within buffer_width of the added frontier cells
    for cell in frontier_cells_added:
        for buffer_cell in project_buffer(cell, buffer_width, rows, cols):
            tracked_cells[buffer_cell] = tracked_cells.get(buffer_cell, 0) + 1
    # Decrement reference counters for all cells within buffer_width of the dropped frontier cells
    for cell in frontier_cells_dropped:
        for buffer_cell in project_buffer(cell, buffer_width, rows, cols):
            tracked_cells[buffer_cell] -= 1
            if tracked_cells[buffer_cell] == 0:
                tracked_cells.pop(buffer_cell)
    # Return updated tracked cells
    return tracked_cells
# phi-field-perimeter-tracking ends here
# [[file:../../org/pyretechnics.org::spread-phi-field][spread-phi-field]]
import numpy as np
import pyretechnics.conversion as conv
from pyretechnics.spot_fire import expected_firebrand_production, spread_firebrands
import pyretechnics.vector_utils as vu


# TODO: Move to pyretechnics.conversion
fire_type_codes = {
    "unburned"      : 0,
    "surface"       : 1,
    "passive_crown" : 2,
    "active_crown"  : 3,
}


# TODO: Move to pyretechnics.vector_utils
def spread_direction_vector_to_angle(vector_3d):
    """
    TODO: Add docstring
    """
    (x, y)       = vu.to_horizontal_plane(vector_3d)
    (r, azimuth) = conv.cartesian_to_azimuthal(x, y)
    return azimuth


def spread_fire_one_timestep(space_time_cubes, output_matrices, frontier_cells, tracked_cells, cube_resolution,
                             start_time, max_timestep, use_wind_limit=True, surface_lw_ratio_model="rothermel",
                             crown_max_lw_ratio=None, max_cells_per_timestep=0.4, buffer_width=3,
                             spot_ignitions={}, spot_config=None, random_generator=None):
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
    # Extract simulation dimensions
    (band_duration, cell_height, cell_width) = cube_resolution

    # Unpack space_time_cubes
    slope_cube  = space_time_cubes["slope"]
    aspect_cube = space_time_cubes["aspect"]

   # Unpack output_matrices
    phi_matrix                = output_matrices["phi"]
    phi_star_matrix           = output_matrices["phi_star"]
    fire_type_matrix          = output_matrices["fire_type"]
    spread_rate_matrix        = output_matrices["spread_rate"]
    spread_direction_matrix   = output_matrices["spread_direction"]
    fireline_intensity_matrix = output_matrices["fireline_intensity"]
    flame_length_matrix       = output_matrices["flame_length"]
    time_of_arrival_matrix    = output_matrices["time_of_arrival"]

    # Extract the firebrands_per_unit_heat from the spot_config if provided
    firebrands_per_unit_heat = spot_config["firebrands_per_unit_heat"] if spot_config else None

    # Initialize max spread rates in the x and y dimensions to 0.0
    max_spread_rate_x = 0.0
    max_spread_rate_y = 0.0

    # Create an empty dictionary to store intermediate fire behavior values per cell
    fire_behavior_dict = {}

    # Compute fire behavior values at start_time and identify the max spread rates in the x and y dimensions
    t0 = int(start_time // band_duration)
    for cell_index in tracked_cells:
        # Unpack cell_index
        (y, x) = cell_index

        # Calculate phi gradient on the horizontal plane
        phi_gradient_xy  = calc_phi_gradient_approx(phi_matrix, cell_width, cell_height, x, y)
        phi_magnitude_xy = vu.vector_magnitude(phi_gradient_xy)

        # Calculate the fire behavior normal to the fire front on the slope-tangential plane
        fire_behavior = burn_cell_toward_phi_gradient(space_time_cubes, (t0, y, x),
                                                      phi_gradient_xy, use_wind_limit,
                                                      surface_lw_ratio_model,
                                                      crown_max_lw_ratio)

        # Check whether cell has a positive phi magnitude
        if phi_magnitude_xy > 0.0:
            # Keep a running tally of the max horizontal spread rates in the x and y dimensions
            (dphi_dx, dphi_dy) = phi_gradient_xy
            phi_magnitude_xy_2 = phi_magnitude_xy ** 2.0
            dphi_dt            = fire_behavior["dphi_dt"]
            spread_rate_x      = -dphi_dt * dphi_dx / phi_magnitude_xy_2
            spread_rate_y      = -dphi_dt * dphi_dy / phi_magnitude_xy_2
            max_spread_rate_x  = max(max_spread_rate_x, abs(spread_rate_x))
            max_spread_rate_y  = max(max_spread_rate_y, abs(spread_rate_y))

            # Integrate the Superbee flux limited phi gradient to make dphi_dt numerically stable
            phi_gradient_xy_limited = calc_phi_gradient(phi_matrix, dphi_dx, dphi_dy, cell_width, cell_height, x, y)
            fire_behavior["dphi_dt"] *= np.dot(phi_gradient_xy, phi_gradient_xy_limited) / phi_magnitude_xy_2

        # Store fire behavior values for later use
        fire_behavior_dict[cell_index] = fire_behavior

    # Calculate timestep using the CFL condition
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
    stop_time = start_time + dt

    # Update the tracked cell values in phi_star_matrix
    for cell_index in tracked_cells:
        (y, x)  = cell_index
        dphi_dt = fire_behavior_dict[cell_index]["dphi_dt"]
        if dphi_dt != 0.0:
            phi_star_matrix[y,x] += dphi_dt * dt

    # Compute fire behavior values at stop_time and update the output_matrices
    t1 = int(stop_time // band_duration)
    for cell_index in tracked_cells:
        # Unpack cell_index
        (y, x) = cell_index

        # Calculate phi gradient on the horizontal plane
        phi_gradient_xy_star  = calc_phi_gradient_approx(phi_star_matrix, cell_width, cell_height, x, y)
        phi_magnitude_xy_star = vu.vector_magnitude(phi_gradient_xy_star)

        # Calculate the fire behavior normal to the fire front on the slope-tangential plane
        fire_behavior_star = burn_cell_toward_phi_gradient(space_time_cubes, (t1, y, x),
                                                           phi_gradient_xy_star, use_wind_limit,
                                                           surface_lw_ratio_model,
                                                           crown_max_lw_ratio)

        # Check whether cell has a positive phi magnitude
        if phi_magnitude_xy_star > 0.0:
            # Integrate the Superbee flux limited phi gradient to make dphi_dt numerically stable
            phi_gradient_xy_star_limited = calc_phi_gradient(phi_star_matrix, *phi_gradient_xy_star,
                                                             cell_width, cell_height, x, y)
            fire_behavior_star["dphi_dt"] *= (np.dot(phi_gradient_xy_star, phi_gradient_xy_star_limited) /
                                              (phi_magnitude_xy_star ** 2.0))

        # Calculate the new phi value at stop_time as phi_next
        fire_behavior     = fire_behavior_dict[cell_index]
        dphi_dt_estimate1 = fire_behavior["dphi_dt"]
        dphi_dt_estimate2 = fire_behavior_star["dphi_dt"]
        dphi_dt_average   = (dphi_dt_estimate1 + dphi_dt_estimate2) / 2.0
        if dphi_dt_average != 0.0:
            phi      = phi_matrix[y,x]
            phi_next = phi + dphi_dt_average * dt

            # Update the tracked cell values in phi_matrix
            phi_matrix[y,x] = phi_next

            # Record fire behavior values in the output_matrices for cells that are burned in this timestep
            # NOTE: This records the fire behavior values at start_time and not at the time of arrival.
            if phi > 0.0 and phi_next <= 0.0:
                fire_type_matrix[y,x]          = fire_type_codes[fire_behavior["fire_type"]]
                spread_rate_matrix[y,x]        = fire_behavior["spread_rate"]
                spread_direction_matrix[y,x]   = spread_direction_vector_to_angle(fire_behavior["spread_direction"])
                fireline_intensity_matrix[y,x] = fire_behavior["fireline_intensity"]
                flame_length_matrix[y,x]       = fire_behavior["flame_length"]
                time_of_arrival_matrix[y,x]    = start_time + dt * phi / (phi - phi_next)

                # Cast firebrands, update firebrand_count_matrix, and update spot_ignitions
                if spot_config:
                    t_cast                   = int(time_of_arrival_matrix[y,x] // band_duration)
                    space_time_coordinate    = (t_cast, y, x)
                    slope                    = slope_cube.get(t_cast, y, x)
                    aspect                   = aspect_cube.get(t_cast, y, x)
                    elevation_gradient       = calc_elevation_gradient(slope, aspect)
                    expected_firebrand_count = expected_firebrand_production(fire_behavior,
                                                                             elevation_gradient,
                                                                             cube_resolution,
                                                                             firebrands_per_unit_heat)
                    new_ignitions            = spread_firebrands(space_time_cubes, output_matrices, cube_resolution,
                                                                 space_time_coordinate, random_generator,
                                                                 expected_firebrand_count, spot_config)
                    if new_ignitions:
                        (ignition_time, ignited_cells) = new_ignitions
                        concurrent_ignited_cells       = spot_ignitions.get(ignition_time)
                        if concurrent_ignited_cells:
                            spot_ignitions[ignition_time] = set.union(ignited_cells, concurrent_ignited_cells)
                        else:
                            spot_ignitions[ignition_time] = ignited_cells

    # Update phi_matrix and time_of_arrival matrix for all cells that ignite a new spot fire before stop_time
    for ignition_time in sorted(spot_ignitions):
        if ignition_time < stop_time:
            ignited_cells = spot_ignitions.pop(ignition_time)
            for cell_index in ignited_cells:
                (y, x) = cell_index
                if phi_matrix[y,x] > 0.0:
                    phi_matrix[y,x]             = -1.0
                    time_of_arrival_matrix[y,x] = ignition_time # FIXME: REVIEW Should I use stop_time instead?
                    tracked_cells[cell_index]   = tracked_cells.get(cell_index, 0)
                    # FIXME: I need to calculate and store the fire_behavior values for these cells

    # Save the new phi_matrix values in phi_star_matrix
    for (y,x) in tracked_cells:
        phi_star_matrix[y,x] = phi_matrix[y,x]

    # Update the sets of frontier cells and tracked cells based on the updated phi matrix
    frontier_cells_new = identify_frontier_cells(phi_matrix, tracked_cells)
    tracked_cells_new  = update_tracked_cells(tracked_cells, frontier_cells, frontier_cells_new,
                                              buffer_width, phi_matrix.shape)

    # Return the updated world state
    return {
        "simulation_time" : stop_time,
        "output_matrices" : output_matrices,
        "frontier_cells"  : frontier_cells_new,
        "tracked_cells"   : tracked_cells_new,
        "spot_ignitions"  : spot_ignitions,
        "random_generator": random_generator,
    }


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
    frontier_cells = identify_frontier_cells(phi_matrix)
    tracked_cells  = identify_tracked_cells(frontier_cells, buffer_width, (rows, cols))

    # Make a copy of the phi matrix to use for intermediate calculations in each timestep
    output_matrices["phi_star"] = np.copy(phi_matrix)

    # Create a numpy.random.Generator object to produce random samples if spot_config is provided
    random_generator = np.random.default_rng(seed=spot_config["random_seed"]) if spot_config else None

    # Spread the fire until an exit condition is reached
    # FIXME: I don't think the "no burnable cells" condition can ever be met currently.
    simulation_time = start_time
    while(simulation_time < max_stop_time and (len(tracked_cells) > 0 or len(spot_ignitions) > 0)):
        # Compute max_timestep based on the remaining time in the temporal band and simulation
        remaining_time_in_band       = band_duration - simulation_time % band_duration
        remaining_time_in_simulation = max_stop_time - simulation_time
        max_timestep                 = min(remaining_time_in_band, remaining_time_in_simulation)

        # Spread fire one timestep
        results = spread_fire_one_timestep(space_time_cubes, output_matrices, frontier_cells, tracked_cells,
                                           cube_resolution, simulation_time, max_timestep, use_wind_limit,
                                           surface_lw_ratio_model, crown_max_lw_ratio, max_cells_per_timestep,
                                           buffer_width, spot_ignitions, spot_config, random_generator)

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
    return {
        "stop_time"      : simulation_time,
        "stop_condition" : "max duration reached" if len(tracked_cells) > 0 else "no burnable cells",
        "output_matrices": output_matrices,
    } | ({
        "spot_ignitions"  : spot_ignitions,
        "random_generator": random_generator,
    } if spot_config else {})
# spread-phi-field ends here
