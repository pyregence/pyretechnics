# [[file:../../org/pyretechnics.org::eulerian-level-set-imports][eulerian-level-set-imports]]
import cython
import cython as cy
import numpy as np
from sortedcontainers import SortedDict
if cython.compiled:
    from cython.cimports.numpy import ndarray
    from cython.cimports.libc.stdlib import malloc, realloc, free
    from cython.cimports.libc.math import pi, floor, sqrt, pow, atan
    from cython.cimports.pyretechnics.cy_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, \
        fclaarr, FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior, SpotConfig, \
        PartialedEllWavelet, CellInputs, EllipticalInfo, Pass1CellOutput
    from cython.cimports.pyretechnics.random import BufferedRandGen
    from cython.cimports.pyretechnics.space_time_cube import ISpaceTimeCube, SpaceTimeCube
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.vector_utils as vu
    import cython.cimports.pyretechnics.fuel_models as fm
    import cython.cimports.pyretechnics.surface_fire as sf
    import cython.cimports.pyretechnics.crown_fire as cf
    import cython.cimports.pyretechnics.spot_fire as spot
    import cython.cimports.pyretechnics.narrow_band_tracking as nbt
else:
    # TODO: Create equivalent Python functions for malloc, realloc, free
    from numpy import ndarray
    from math import pi, floor, sqrt, pow, atan
    from pyretechnics.py_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, \
        fclaarr, FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior, SpotConfig, \
        PartialedEllWavelet, CellInputs, EllipticalInfo, Pass1CellOutput
    from pyretechnics.random import BufferedRandGen
    from pyretechnics.space_time_cube import ISpaceTimeCube, SpaceTimeCube
    import pyretechnics.conversion as conv
    import pyretechnics.vector_utils as vu
    import pyretechnics.fuel_models as fm
    import pyretechnics.surface_fire as sf
    import pyretechnics.crown_fire as cf
    import pyretechnics.spot_fire as spot
    import pyretechnics.narrow_band_tracking as nbt
# eulerian-level-set-imports ends here
# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients-approx][phi-field-spatial-gradients-approx]]
# NOTE: No longer used in tight loops.
@cy.cfunc
@cy.exceptval(check=False)
def calc_dphi_dx_approx(phi_matrix: cy.float[:,::1], dx: cy.float, x: pyidx, y: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of `phi_matrix` in the x (west->east)
    direction at grid cell (x,y) given the cell width `dx`.
    """
    east_x: pyidx = x + 1
    west_x: pyidx = x - 1
    return (phi_matrix[2+y, 2+east_x] - phi_matrix[2+y, 2+west_x]) / (2.0 * dx)


# NOTE: No longer used in tight loops.
@cy.cfunc
@cy.exceptval(check=False)
def calc_dphi_dy_approx(phi_matrix: cy.float[:,::1], dy: cy.float, x: pyidx, y: pyidx) -> cy.float:
    """
    Calculate the spatial gradient of `phi_matrix` in the y (south->north)
    direction at grid cell (x,y) given the cell height `dy`.
    """
    north_y: pyidx = y + 1
    south_y: pyidx = y - 1
    return (phi_matrix[2+north_y, 2+x] - phi_matrix[2+south_y, 2+x]) / (2.0 * dy)


# NOTE: No longer used in tight loops.
@cy.cfunc
@cy.exceptval(check=False)
def calc_phi_gradient_approx(phi_matrix: cy.float[:,::1], dx: cy.float, dy: cy.float, x: pyidx, y: pyidx) -> vec_xy:
    """
    Calculate the spatial gradient of `phi_matrix` at grid cell (x,y)
    given the cell width `dx` and the cell height `dy`.
    """
    dphi_dx: cy.float = calc_dphi_dx_approx(phi_matrix, dx, x, y)
    dphi_dy: cy.float = calc_dphi_dy_approx(phi_matrix, dy, x, y)
    return (dphi_dx, dphi_dy)
# phi-field-spatial-gradients-approx ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector][phi-field-normal-vector]]
# TODO: Remove unused function
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_phi_normal_vector(phi_matrix: cy.float[:,::1], dx: cy.float, dy: cy.float, x: pyidx, y: pyidx) -> vec_xy:
    """
    Calculate the phi field normal vector in the x and y dimensions.

    - n_x: eastward component of the unit normal vector
    - n_y: northward component of the unit normal vector
    """
    return vu.as_unit_vector_2d(calc_phi_gradient_approx(phi_matrix, dx, dy, x, y)) # (n_x, n_y)
# phi-field-normal-vector ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector-angle][phi-field-normal-vector-angle]]
# TODO: Remove unused function
@cy.cfunc
@cy.exceptval(check=False)
def calc_phi_normal_azimuth(phi_normal_vector: vec_xy) -> cy.float:
    """
    Calculate the angle (measured in degrees clockwise from North)
    to which the phi field's normal vector points.
    """
    (n_x, n_y) = phi_normal_vector
    angle: cy.float
    if n_x > 0.0:
        if n_y >= 0.0:
            angle = 0.5 * pi - atan(n_y / n_x)
        elif n_y < 0.0:
            angle = 0.5 * pi + atan(abs(n_y) / n_x)
    elif n_x < 0.0:
        if n_y >= 0.0:
            angle = 1.5 * pi + atan(n_y / abs(n_x))
        elif n_y < 0.0:
            angle = 1.5 * pi - atan(n_y / n_x)
    else:
        if n_y >= 0.0:
            angle = 0.0
        elif n_y < 0.0:
            angle = pi
    return conv.rad_to_deg(angle)
# phi-field-normal-vector-angle ends here
# [[file:../../org/pyretechnics.org::superbee-flux-limiter][superbee-flux-limiter]]
@cy.cfunc
@cy.exceptval(check=False)
def half_superbee_dphi_up(dphi_up: cy.float, dphi_loc: cy.float) -> cy.float:
    """
    Logically like calc_superbee_flux_limiter() but returns a result multiplied by (0.5 * dphi_loc).

    NOTE: This is more numerically stable than calc_superbee_flux_limiter().
    """
    s_loc             : cy.float = 1.0 if dphi_loc >= 0.0 else -1.0
    are_opposite_signs: cy.bint  = (s_loc * dphi_up) <= 0.0
    if are_opposite_signs:
        return 0.0
    a_up : cy.float = abs(dphi_up)
    a_loc: cy.float = abs(dphi_loc)
    return s_loc * max(min(a_up / 2.0, a_loc),
                       min(a_up, a_loc / 2.0))
# superbee-flux-limiter ends here
# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients][phi-field-spatial-gradients]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_dphi_flim_x(p00: cy.float, pw2: cy.float, pw1: cy.float, pe1: cy.float, pe2: cy.float) -> cy.float:
    dphi_up : cy.float
    dphi_loc: cy.float
    phi_east: cy.float
    phi_west: cy.float

    dphi_loc = pe1 - p00
    if pe1 >= pw1:
        dphi_up  = p00 - pw1
        phi_east = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up  = pe2 - pe1
        phi_east = pe1 - half_superbee_dphi_up(dphi_up, dphi_loc)

    dphi_loc = pw1 - p00
    if pe1 >= pw1:
        dphi_up  = pw2 - pw1
        phi_west = pw1 - half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up  = p00 - pe1
        phi_west = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)
    return (phi_east - phi_west)


# NOTE: This is actually the same function as the previous one. But
#       who knows, maybe we get a performance gain by differentiating
#       code sites.
@cy.cfunc
@cy.exceptval(check=False)
def calc_dphi_flim_y(p00: cy.float, ps2: cy.float, ps1: cy.float, pn1: cy.float, pn2: cy.float) -> cy.float:
    dphi_up  : cy.float
    dphi_loc : cy.float
    phi_north: cy.float
    phi_south: cy.float

    dphi_loc = pn1 - p00
    if pn1 >= ps1:
        dphi_up   = p00 - ps1
        phi_north = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up   = pn2 - pn1
        phi_north = pn1 - half_superbee_dphi_up(dphi_up, dphi_loc)

    dphi_loc = ps1 - p00
    if pn1 >= ps1:
        dphi_up   = ps2 - ps1
        phi_south = ps1 - half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up   = p00 - pn1
        phi_south = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)

    return (phi_north - phi_south)
# phi-field-spatial-gradients ends here
# [[file:../../org/pyretechnics.org::calc-fireline-normal-behavior][calc-fireline-normal-behavior]]
# TODO: Move these to a shared module and use throughout the literate program
# NOTE: It would be better to use a cython enum here, but that's not supported in pure python syntax.
fire_type_unburned      = cy.declare(cy.int, 0)
fire_type_surface       = cy.declare(cy.int, 1)
fire_type_crown_passive = cy.declare(cy.int, 2)
fire_type_crown_active  = cy.declare(cy.int, 3)


# TODO: Move this to pyretechnics.vector_utils and use throughout the literate program
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def calc_elevation_gradient(slope: cy.float, aspect: cy.float) -> vec_xy:
    """
    Returns the elevation gradient (dz_dx: rise/run, dz_dy: rise/run) given:
    - slope  :: rise/run
    - aspect :: degrees clockwise from North
    """
    return conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))


@cy.cfunc
@cy.exceptval(check=False)
def calc_phi_gradient_on_slope(phi_gradient_xy: vec_xy, elevation_gradient: vec_xy) -> vec_xyz:
    """
    Returns the gradient of phi projected onto the slope-tangential plane as a 3D (x,y,z) vector (in phi/m) given:
    - phi_gradient_xy    :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - elevation_gradient :: (dz_dx: rise/run, dz_dy: rise/run)
    """
    (dphi_dx, dphi_dy)        = phi_gradient_xy
    phi_gradient_xyz: vec_xyz = (dphi_dx, dphi_dy, 0.0)
    if vu.vector_magnitude_2d(elevation_gradient) == 0.0:
        return phi_gradient_xyz
    else:
        slope_normal_vector: vec_xyz  = vu.get_slope_normal_vector(elevation_gradient) # (x,y,z) unit vector
        phi_slope_agreement: cy.float = vu.dot_3d(phi_gradient_xyz, slope_normal_vector)
        dphi_dx_on_slope   : cy.float = phi_gradient_xyz[0] - phi_slope_agreement * slope_normal_vector[0]
        dphi_dy_on_slope   : cy.float = phi_gradient_xyz[1] - phi_slope_agreement * slope_normal_vector[1]
        dphi_dz_on_slope   : cy.float = phi_gradient_xyz[2] - phi_slope_agreement * slope_normal_vector[2]
        return (dphi_dx_on_slope, dphi_dy_on_slope, dphi_dz_on_slope)


# FIXME: Do I switch to cruz_passive_crown_fire_spread_rate() if the normal_spread_rate < critical_spread_rate?
#        Did I do this correctly in calc_crown_fire_behavior_in_direction?
@cy.cfunc
@cy.exceptval(check=False)
def calc_fireline_normal_behavior(fire_behavior_max: FireBehaviorMax, phi_gradient_xyz: vec_xyz) -> SpreadBehavior:
    """
    Given these inputs:
    - fire_behavior_max  :: FireBehaviorMax struct of max surface or crown fire behavior values
      - max_fire_type          :: 0 (unburned), 1 (surface), 2 (passive_crown), or 3 (active_crown)
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_fireline_intensity :: kW/m
      - max_flame_length       :: m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
      - critical_spread_rate   :: m/min (Required for crown fires only)
    - phi_gradient_xyz   :: (dphi_dx: phi/m, dphi_dy: phi/m, dphi_dz: phi/m) 3D vector on the slope-tangential plane

    return a SpreadBehavior struct containing these keys:
    - dphi_dt            :: phi/min (on the slope-tangential plane)
    - fire_type          :: 0 (unburned), 1 (surface), 2 (passive_crown), or 3 (active_crown)
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m

    NOTE: This function should work for surface or crown fires interchangeably.
    """
    #================================================================================================
    # Calculate the magnitude of the phi gradient
    #================================================================================================

    phi_magnitude: cy.float = vu.vector_magnitude_3d(phi_gradient_xyz) # phi/m

    #================================================================================================
    # Check whether cell is on the fire perimeter and burning
    #================================================================================================

    if (phi_magnitude == 0.0 or fire_behavior_max.max_spread_rate == 0.0):
        # This location is not on the fire perimeter and/or is not burning

        #================================================================================================
        # Set the spread direction to the phi gradient direction, max spread direction, upslope, or North
        #================================================================================================

        spread_direction: vec_xyz = (vu.scale_3d(1.0 / phi_magnitude, phi_gradient_xyz)
                                     if phi_magnitude > 0.0
                                     else fire_behavior_max.max_spread_direction)

        #============================================================================================
        # Return zero surface/crown fire behavior
        #============================================================================================

        return SpreadBehavior(
            dphi_dt            = 0.0,
            fire_type          = fire_type_unburned,
            spread_rate        = 0.0,
            spread_direction   = spread_direction,
            fireline_intensity = 0.0,
            flame_length       = 0.0,
        )

    else:
        # This location is on the fire perimeter and is burning

        #============================================================================================
        # Unpack the fire_behavior_max struct
        #============================================================================================

        heading_fire_type         : cy.int   = fire_behavior_max.max_fire_type
        heading_spread_rate       : cy.float = fire_behavior_max.max_spread_rate        # m/min
        heading_spread_direction  : vec_xyz  = fire_behavior_max.max_spread_direction   # (x,y,z) unit vector
        heading_spread_vector     : vec_xyz  = vu.scale_3d(heading_spread_rate,
                                                           heading_spread_direction)    # (x,y,z) m/min vector
        heading_fireline_intensity: cy.float = fire_behavior_max.max_fireline_intensity # kW/m
        length_to_width_ratio     : cy.float = fire_behavior_max.length_to_width_ratio  # unitless
        eccentricity              : cy.float = fire_behavior_max.eccentricity           # unitless
        critical_spread_rate      : cy.float = fire_behavior_max.critical_spread_rate   # m/min

        #============================================================================================
        # Calculate the backing and flanking fire spread rates
        #============================================================================================

        backing_adjustment  : cy.float = (1.0 - eccentricity) / (1.0 + eccentricity)  # unitless
        backing_spread_rate : cy.float = heading_spread_rate * backing_adjustment     # m/min
        flanking_spread_rate: cy.float = ((heading_spread_rate + backing_spread_rate)
                                          / (2.0 * length_to_width_ratio))            # m/min

        #============================================================================================
        # Calculate dphi/dt
        #============================================================================================

        A      : cy.float = (heading_spread_rate - backing_spread_rate) / (2.0 * heading_spread_rate) # unitless
        B      : cy.float = vu.dot_3d(heading_spread_vector, phi_gradient_xyz)                        # phi/min
        C      : cy.float = flanking_spread_rate / heading_spread_rate                                # unitless
        D      : cy.float = pow((heading_spread_rate * phi_magnitude), 2.0)                           # (phi/min)^2
        E      : cy.float = ((length_to_width_ratio * length_to_width_ratio) - 1.0) * (B * B)         # (phi/min)^2
        dphi_dt: cy.float = -(A * B + C * sqrt(D + E))                                                # phi/min

        #============================================================================================
        # Calculate fire behavior normal to the fire perimeter
        #============================================================================================

        normal_spread_rate       : cy.float = -dphi_dt / phi_magnitude                        # m/min
        normal_spread_direction  : vec_xyz  = vu.as_unit_vector_3d(phi_gradient_xyz)          # (x,y,z) unit vector
        normal_adjustment        : cy.float = normal_spread_rate / heading_spread_rate        # unitless
        normal_fireline_intensity: cy.float = heading_fireline_intensity * normal_adjustment  # kW/m
        normal_flame_length      : cy.float = sf.calc_flame_length(normal_fireline_intensity) # m
        normal_fire_type         : cy.int   = (fire_type_surface if heading_fire_type == fire_type_surface
                                               else fire_type_crown_active if normal_spread_rate > critical_spread_rate
                                               else fire_type_crown_passive)

        #========================================================================================
        # Return the surface/crown fire behavior normal to the fire perimeter
        #========================================================================================

        return SpreadBehavior(
            dphi_dt            = dphi_dt,                   # phi/min
            fire_type          = normal_fire_type,          # surface, passive_crown, or active_crown
            spread_rate        = normal_spread_rate,        # m/min
            spread_direction   = normal_spread_direction,   # (x,y,z) unit vector
            fireline_intensity = normal_fireline_intensity, # kW/m
            flame_length       = normal_flame_length,       # m
        )
# calc-fireline-normal-behavior ends here
# [[file:../../org/pyretechnics.org::burn-cell-toward-phi-gradient][burn-cell-toward-phi-gradient]]
# TODO: Turn this into a struct once its methods have been removed
@cy.cclass
class SpreadInputs:
    """
    A fast-access data structure for reading inputs in performance-critical code.
    """
    rows                         : pyidx
    cols                         : pyidx
    band_duration                : cy.float # minutes
    cell_height                  : cy.float # meters
    cell_width                   : cy.float # meters
    slope                        : ISpaceTimeCube
    aspect                       : ISpaceTimeCube
    fuel_model                   : ISpaceTimeCube
    canopy_cover                 : ISpaceTimeCube
    canopy_height                : ISpaceTimeCube
    canopy_base_height           : ISpaceTimeCube
    canopy_bulk_density          : ISpaceTimeCube
    wind_speed_10m               : ISpaceTimeCube
    upwind_direction             : ISpaceTimeCube
    fuel_moisture_dead_1hr       : ISpaceTimeCube
    fuel_moisture_dead_10hr      : ISpaceTimeCube
    fuel_moisture_dead_100hr     : ISpaceTimeCube
    fuel_moisture_live_herbaceous: ISpaceTimeCube
    fuel_moisture_live_woody     : ISpaceTimeCube
    foliar_moisture              : ISpaceTimeCube
    temperature                  : ISpaceTimeCube
    fuel_spread_adjustment       : ISpaceTimeCube
    weather_spread_adjustment    : ISpaceTimeCube
    fuel_model_cache             : cy.pointer(FuelModel)


    def __cinit__(self,
                  cube_resolution              : tuple[cy.float, cy.float, cy.float],
                  slope                        : ISpaceTimeCube,
                  aspect                       : ISpaceTimeCube,
                  fuel_model                   : ISpaceTimeCube,
                  canopy_cover                 : ISpaceTimeCube,
                  canopy_height                : ISpaceTimeCube,
                  canopy_base_height           : ISpaceTimeCube,
                  canopy_bulk_density          : ISpaceTimeCube,
                  wind_speed_10m               : ISpaceTimeCube,
                  upwind_direction             : ISpaceTimeCube,
                  fuel_moisture_dead_1hr       : ISpaceTimeCube,
                  fuel_moisture_dead_10hr      : ISpaceTimeCube,
                  fuel_moisture_dead_100hr     : ISpaceTimeCube,
                  fuel_moisture_live_herbaceous: ISpaceTimeCube,
                  fuel_moisture_live_woody     : ISpaceTimeCube,
                  foliar_moisture              : ISpaceTimeCube,
                  temperature                  : ISpaceTimeCube,
                  fuel_spread_adjustment       : ISpaceTimeCube,
                  weather_spread_adjustment    : ISpaceTimeCube) -> cy.void:
        (_bands, rows, cols)               = slope.shape
        self.rows                          = rows
        self.cols                          = cols
        self.band_duration                 = cube_resolution[0]
        self.cell_height                   = cube_resolution[1]
        self.cell_width                    = cube_resolution[2]
        self.slope                         = slope
        self.aspect                        = aspect
        self.fuel_model                    = fuel_model
        self.canopy_cover                  = canopy_cover
        self.canopy_height                 = canopy_height
        self.canopy_base_height            = canopy_base_height
        self.canopy_bulk_density           = canopy_bulk_density
        self.wind_speed_10m                = wind_speed_10m
        self.upwind_direction              = upwind_direction
        self.fuel_moisture_dead_1hr        = fuel_moisture_dead_1hr
        self.fuel_moisture_dead_10hr       = fuel_moisture_dead_10hr
        self.fuel_moisture_dead_100hr      = fuel_moisture_dead_100hr
        self.fuel_moisture_live_herbaceous = fuel_moisture_live_herbaceous
        self.fuel_moisture_live_woody      = fuel_moisture_live_woody
        self.foliar_moisture               = foliar_moisture
        self.temperature                   = temperature
        self.fuel_spread_adjustment        = fuel_spread_adjustment
        self.weather_spread_adjustment     = weather_spread_adjustment
        self.__init_fuel_models()


    # TODO: Move this code to fuel_models.py
    @cy.cfunc
    def __init_fuel_models(self) -> cy.void:
        # Allocate an empty FuelModel array in memory
        fuel_model_cache: cy.pointer(FuelModel) = cy.cast(cy.pointer(FuelModel),
                                                          malloc(300 * cy.sizeof(FuelModel)))
        # Verify that it was created
        if not fuel_model_cache:
            # Something went wrong with malloc
            raise MemoryError()
        else:
            # Copy FuelModels from fm.fuel_model_table into fuel_model_cache
            fuel_model: FuelModel
            for fuel_model in fm.fuel_model_table.values():
                fuel_model_cache[fuel_model.number] = fuel_model
            # Save fuel_model_cache in the SpreadInputs object
            self.fuel_model_cache = fuel_model_cache


    # TODO: Inline this code at its call sites
    @cy.cfunc
    @cy.inline
    @cy.exceptval(check=False)
    def get_fm_struct(self, fm_number: pyidx) -> FuelModel:
        return self.fuel_model_cache[fm_number]


    def __dealloc__(self) -> cy.void:
        free(self.fuel_model_cache) # no-op if self.fuel_model_cache is NULL


# TODO: Inline this code at its call sites
@cy.cfunc
@cy.inline
def make_SpreadInputs(cube_shape      : tuple[pyidx, pyidx, pyidx],
                      cube_resolution : tuple[cy.float, cy.float, cy.float],
                      space_time_cubes: dict) -> SpreadInputs:
    return SpreadInputs(cube_resolution,
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
                        space_time_cubes.get("temperature", SpaceTimeCube(cube_shape, -1.0)),
                        space_time_cubes.get("fuel_spread_adjustment", SpaceTimeCube(cube_shape, 1.0)),
                        space_time_cubes.get("weather_spread_adjustment", SpaceTimeCube(cube_shape, 1.0)))


@cy.cfunc
@cy.exceptval(check=False)
def lookup_cell_inputs(spread_inputs: SpreadInputs, space_time_coordinate: coord_tyx) -> CellInputs:
    """
    Reads the inputs for a given cell from the space-time cubes, returning a `CellInputs` struct.
    """
    # Unpack the space_time_coordinate
    t: pyidx = space_time_coordinate[0]
    y: pyidx = space_time_coordinate[1]
    x: pyidx = space_time_coordinate[2]

    # Topography, Fuel Model, and Vegetation
    slope              : cy.float = spread_inputs.slope.get(t, y, x)               # rise/run
    aspect             : cy.float = spread_inputs.aspect.get(t, y, x)              # degrees clockwise from North
    fuel_model_number  : cy.float = spread_inputs.fuel_model.get(t, y, x)          # integer index in fm.fuel_model_table
    canopy_cover       : cy.float = spread_inputs.canopy_cover.get(t, y, x)        # 0-1
    canopy_height      : cy.float = spread_inputs.canopy_height.get(t, y, x)       # m
    canopy_base_height : cy.float = spread_inputs.canopy_base_height.get(t, y, x)  # m
    canopy_bulk_density: cy.float = spread_inputs.canopy_bulk_density.get(t, y, x) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m               : cy.float = spread_inputs.wind_speed_10m.get(t, y, x)                # km/hr
    upwind_direction             : cy.float = spread_inputs.upwind_direction.get(t, y, x)              # degrees clockwise from North
    fuel_moisture_dead_1hr       : cy.float = spread_inputs.fuel_moisture_dead_1hr.get(t, y, x)        # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr      : cy.float = spread_inputs.fuel_moisture_dead_10hr.get(t, y, x)       # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr     : cy.float = spread_inputs.fuel_moisture_dead_100hr.get(t, y, x)      # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous: cy.float = spread_inputs.fuel_moisture_live_herbaceous.get(t, y, x) # kg moisture/kg ovendry weight
    fuel_moisture_live_woody     : cy.float = spread_inputs.fuel_moisture_live_woody.get(t, y, x)      # kg moisture/kg ovendry weight
    foliar_moisture              : cy.float = spread_inputs.foliar_moisture.get(t, y, x)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment   : cy.float = (spread_inputs.fuel_spread_adjustment.get(t, y, x)
                                           if spread_inputs.fuel_spread_adjustment is not None
                                           else 1.0)                                         # float >= 0.0
    weather_spread_adjustment: cy.float = (spread_inputs.weather_spread_adjustment.get(t, y, x)
                                           if spread_inputs.weather_spread_adjustment is not None
                                           else 1.0)                                         # float >= 0.0

    # Pack values into a CellInputs struct and return
    return CellInputs(
        slope                         = slope,
        aspect                        = aspect,
        fuel_model_number             = fuel_model_number,
        canopy_cover                  = canopy_cover,
        canopy_height                 = canopy_height,
        canopy_base_height            = canopy_base_height,
        canopy_bulk_density           = canopy_bulk_density,
        wind_speed_10m                = wind_speed_10m,
        upwind_direction              = upwind_direction,
        fuel_moisture_dead_1hr        = fuel_moisture_dead_1hr,
        fuel_moisture_dead_10hr       = fuel_moisture_dead_10hr,
        fuel_moisture_dead_100hr      = fuel_moisture_dead_100hr,
        fuel_moisture_live_herbaceous = fuel_moisture_live_herbaceous,
        fuel_moisture_live_woody      = fuel_moisture_live_woody,
        foliar_moisture               = foliar_moisture,
        fuel_spread_adjustment        = fuel_spread_adjustment,
        weather_spread_adjustment     = weather_spread_adjustment,
    )


# NOTE: We're only going through these annoying calculations because
# we are required to return a spread_direction unit vector, which is
# of questionable value.
#
# IMPROVEMENT: We wouldn't have to go through this trouble if we simply
# returned a Cartesian speed vector instead, which would play more
# nicely with the rest of the code.
@cy.cfunc
@cy.exceptval(check=False)
def unburned_SpreadBehavior(elevation_gradient: vec_xy, phi_gradient_xyz: vec_xyz) -> SpreadBehavior:
    # Create a 3D unit vector pointing...
    spread_direction: vec_xyz
    if (phi_gradient_xyz[0] != 0.0 and phi_gradient_xyz[1] != 0.0 and phi_gradient_xyz[2] != 0.0):
        # ...in the direction of the spatial phi gradient
        spread_direction = vu.as_unit_vector_3d(phi_gradient_xyz)
    elif (elevation_gradient[0] != 0.0 and elevation_gradient[1] != 0.0):
        # ...upslope on the slope-tangential plane
        slope_vector_3d: vec_xyz = vu.to_slope_plane(elevation_gradient, elevation_gradient)
        spread_direction         = vu.as_unit_vector_3d(slope_vector_3d)
    else:
        # ...to the North
        spread_direction = (0.0, 1.0, 0.0)

    # Return zero surface fire behavior
    return SpreadBehavior(
        dphi_dt            = 0.0,
        fire_type          = fire_type_unburned,
        spread_rate        = 0.0,
        spread_direction   = spread_direction,
        fireline_intensity = 0.0,
        flame_length       = 0.0,
    )


@cy.cclass
class SpreadState:
    """
    Stores the stateful data associated with a fire spread simulation.
    """
    cube_shape        : tuple[pyidx, pyidx, pyidx]
    phi               : cy.float[:,::1] # 2D float array of values in [-1,1]
    phi_star          : cy.float[:,::1] # 2D float array of values in [-1,1]
    fire_type         : cy.uchar[:,::1] # 2D byte array (0-3)
    spread_rate       : cy.float[:,::1] # 2D float array (m/min)
    spread_direction  : cy.float[:,::1] # 2D float array (degrees clockwise from North)
    fireline_intensity: cy.float[:,::1] # 2D float array (kW/m)
    flame_length      : cy.float[:,::1] # 2D float array (m)
    time_of_arrival   : cy.float[:,::1] # 2D float array (min)


    # TODO: Initialize output matrices to NaN if possible
    def __init__(self, cube_shape: tuple[pyidx, pyidx, pyidx]) -> cy.void:
        # Extract the grid_shape from the cube_shape
        grid_rows : pyidx               = cube_shape[1]
        grid_cols : pyidx               = cube_shape[2]
        grid_shape: tuple[pyidx, pyidx] = (grid_rows, grid_cols)
        # Create the initial 2D arrays
        # NOTE: The phi matrix is padded by 2 cells on each side to avoid the cost of checking bounds.
        self.cube_shape         = cube_shape
        self.phi                = np.ones((grid_rows + 4, grid_cols + 4), dtype="float32")
        self.phi_star           = np.ones((grid_rows + 4, grid_cols + 4), dtype="float32")
        self.fire_type          = np.zeros(grid_shape, dtype="uint8")
        self.spread_rate        = np.zeros(grid_shape, dtype="float32")
        self.spread_direction   = np.zeros(grid_shape, dtype="float32")
        self.fireline_intensity = np.zeros(grid_shape, dtype="float32")
        self.flame_length       = np.zeros(grid_shape, dtype="float32")
        self.time_of_arrival    = np.full(grid_shape, -1.0, dtype="float32")


    @cy.ccall
    def ignite_cell(self, ignited_cell: coord_yx) -> SpreadState:
        # Extract coords
        y: pyidx = ignited_cell[0]
        x: pyidx = ignited_cell[1]
        # Overwrite phi and phi_star state
        self.phi[2+y,2+x]      = -1.0
        self.phi_star[2+y,2+x] = -1.0
        # Return the updated SpreadState object
        return self


    @cy.ccall
    def ignite_cells(self, lower_left_corner: coord_yx, ignition_matrix: cy.float[:,::1]) -> SpreadState:
        # Extract coords
        rows : pyidx = ignition_matrix.shape[0]
        cols : pyidx = ignition_matrix.shape[1]
        min_y: pyidx = lower_left_corner[0]
        min_x: pyidx = lower_left_corner[1]
        max_y: pyidx = min_y + rows
        max_x: pyidx = min_x + cols
        # Overwrite phi and phi_star state
        self.phi[2+min_y:2+max_y,2+min_x:2+max_x]      = ignition_matrix
        self.phi_star[2+min_y:2+max_y,2+min_x:2+max_x] = ignition_matrix
        # Return the updated SpreadState object
        return self


    @cy.ccall
    def get_burned_matrices(self, layers: list[str]|None = None) -> dict:
        # Find bounding box of burned area
        burned_mask: tuple[ndarray, ndarray] = np.nonzero(self.fire_type)
        min_y      : pyidx                   = burned_mask[0].min()
        max_y      : pyidx                   = burned_mask[0].max()
        min_x      : pyidx                   = burned_mask[1].min()
        max_x      : pyidx                   = burned_mask[1].max()
        # Prepare the 2D arrays in a dict
        available_matrices: dict[str, ndarray] = {
            "phi"               : self.phi,
            "fire_type"         : self.fire_type,
            "spread_rate"       : self.spread_rate,
            "spread_direction"  : self.spread_direction,
            "fireline_intensity": self.fireline_intensity,
            "flame_length"      : self.flame_length,
            "time_of_arrival"   : self.time_of_arrival,
        }
        # Set selected_layers to layers if specified and otherwise to all available layers
        selected_layers: list[str] = layers if layers is not None else list(available_matrices.keys())
        # Clip the 2D arrays from selected_layers to the bounding box
        clipped_matrices: dict[str, ndarray] = {
            k: np.copy(available_matrices[k][min_y:max_y+1, min_x:max_x+1]) for k in selected_layers
        }
        # Return the clipped_matrices along with their lower_left_corner for reference
        return {
            "cube_shape"       : self.cube_shape,
            "lower_left_corner": (min_y, min_x),
            "clipped_matrices" : clipped_matrices,
        }


    @cy.ccall
    def get_full_matrices(self, layers: list[str]|None = None) -> dict:
        # Prepare the 2D arrays in a dict
        available_matrices: dict[str, ndarray] = {
            "phi"               : np.asarray(self.phi),
            "fire_type"         : np.asarray(self.fire_type),
            "spread_rate"       : np.asarray(self.spread_rate),
            "spread_direction"  : np.asarray(self.spread_direction),
            "fireline_intensity": np.asarray(self.fireline_intensity),
            "flame_length"      : np.asarray(self.flame_length),
            "time_of_arrival"   : np.asarray(self.time_of_arrival),
        }
        # Return the matrices in layers if specified and otherwise return all available layers
        if layers is None:
            return available_matrices
        else:
            return {
                k: available_matrices[k] for k in layers
            }


    @cy.ccall
    def copy(self) -> SpreadState:
        # Create an empty SpreadState object
        new_spread_state: SpreadState = SpreadState.__new__(SpreadState, self.cube_shape)
        # Initialize its fields with copies of the base object's fields
        new_spread_state.cube_shape         = self.cube_shape # tuples are immutable
        new_spread_state.phi                = np.copy(self.phi)
        new_spread_state.phi_star           = np.copy(self.phi_star)
        new_spread_state.fire_type          = np.copy(self.fire_type)
        new_spread_state.spread_rate        = np.copy(self.spread_rate)
        new_spread_state.spread_direction   = np.copy(self.spread_direction)
        new_spread_state.fireline_intensity = np.copy(self.fireline_intensity)
        new_spread_state.flame_length       = np.copy(self.flame_length)
        new_spread_state.time_of_arrival    = np.copy(self.time_of_arrival)
        # Return the initialized SpreadState object
        return new_spread_state
# burn-cell-toward-phi-gradient ends here
# [[file:../../org/pyretechnics.org::phi-field-perimeter-tracking][phi-field-perimeter-tracking]]
@cy.cfunc
@cy.inline
def encode_cell_index(y: pyidx, x: pyidx) -> object: # ulong
    """
    Encodes a (y, x) tuple into a single Python integer object.
    This enables a more efficient memory layout than a tuple of Python integers.
    """
    return (cy.cast(cy.ulonglong, y) << 32) + cy.cast(cy.ulonglong, x)


@cy.cfunc
@cy.exceptval(check=False)
def decode_cell_index(encoded_cell_index: object) -> coord_yx:
    cell_index: cy.ulonglong = encoded_cell_index
    y         : pyidx        = cell_index >> 32
    x         : pyidx        = cy.cast(cy.uint, cell_index) # NOTE: faster than (cell_index & 0xFFFFFFFF)
    return (y, x)


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def opposite_phi_signs(phi_matrix: cy.float[:,::1], y1: pyidx, x1: pyidx, y2: pyidx, x2: pyidx) -> cy.bint:
    """
    Return True if the phi values at cells (x1,y1) and (x2,y2) have opposite signs.

    NOTE: phi_matrix has a 2 cell buffer on all 4 sides, so we have to add (2,2) to each cell index.
    """
    return phi_matrix[2+y1, 2+x1] * phi_matrix[2+y2, 2+x2] < 0.0


# TODO: Pass a 2D fuel_model_array instead for speed
@cy.cfunc
@cy.exceptval(check=False)
def is_frontier_cell(phi_matrix     : cy.float[:,::1],
                     fuel_model_cube: ISpaceTimeCube,
                     t              : pyidx,
                     y              : pyidx,
                     x              : pyidx) -> cy.bint:
    # Compare (north, south, east, west) neighboring cell pairs for opposite phi signs
    north_y: pyidx = y+1
    south_y: pyidx = y-1
    east_x : pyidx = x+1
    west_x : pyidx = x-1
    return ((
        # Check north
        opposite_phi_signs(phi_matrix, y, x, north_y, x) and spot.is_burnable_cell(fuel_model_cube, t, north_y, x)
    ) or (
        # Check south
        opposite_phi_signs(phi_matrix, y, x, south_y, x) and spot.is_burnable_cell(fuel_model_cube, t, south_y, x)
    ) or (
        # Check east
        opposite_phi_signs(phi_matrix, y, x, y, east_x) and spot.is_burnable_cell(fuel_model_cube, t, y, east_x)
    ) or (
        # Check west
        opposite_phi_signs(phi_matrix, y, x, y, west_x) and spot.is_burnable_cell(fuel_model_cube, t, y, west_x)
    )) and (
        # Is this cell burnable?
        spot.is_burnable_cell(fuel_model_cube, t, y, x)
    )


# TODO: Is it faster to build up a list or a set?
@cy.cfunc
def identify_all_frontier_cells(phi_matrix     : cy.float[:,::1],
                                fuel_model_cube: ISpaceTimeCube,
                                t              : pyidx,
                                rows           : pyidx,
                                cols           : pyidx) -> set:
    """
    TODO: Add docstring
    """
    frontier_cells: set = set()
    y             : pyidx
    x             : pyidx
    for y in range(rows):
        for x in range(cols):
            if is_frontier_cell(phi_matrix, fuel_model_cube, t, y, x):
                frontier_cells.add(encode_cell_index(y, x))
    return frontier_cells


@cy.cfunc
def identify_tracked_cells(frontier_cells: set,
                           buffer_width  : pyidx,
                           rows          : pyidx,
                           cols          : pyidx) -> nbt.NarrowBandTracker:
    """
    TODO: Add docstring
    """
    tracked_cells     : nbt.NarrowBandTracker = nbt.new_NarrowBandTracker(rows, cols)
    encoded_cell_index: object
    y                 : pyidx
    x                 : pyidx
    for encoded_cell_index in frontier_cells:
        (y, x) = decode_cell_index(encoded_cell_index)
        nbt.inc_square_around(tracked_cells, y, x, buffer_width)
    return tracked_cells
# phi-field-perimeter-tracking ends here
# [[file:../../org/pyretechnics.org::spread-phi-field][spread-phi-field]]
# TODO: OPTIM We might want to pass in the CellInputs and avoid looking up the SpreadInputs again here.
@cy.cfunc
@cy.exceptval(check=False)
def spot_from_burned_cell(spread_inputs   : SpreadInputs,
                          fire_type_matrix: cy.uchar[:,::1],
                          y               : pyidx,
                          x               : pyidx,
                          fire_behavior   : SpreadBehavior,
                          time_of_arrival : cy.float,
                          random_generator: BufferedRandGen,
                          spot_config     : SpotConfig,
                          spot_ignitions  : SortedDict[float, set]) -> cy.void:
    """
    Schedules the future spot ignitions following from burning the given cell.
    Mutates `spot_ignitions` (and `random_generator`).
    """
    t_cast                  : pyidx       = int(time_of_arrival // spread_inputs.band_duration)
    slope                   : cy.float    = spread_inputs.slope.get(t_cast, y, x)
    aspect                  : cy.float    = spread_inputs.aspect.get(t_cast, y, x)
    elevation_gradient      : vec_xy      = calc_elevation_gradient(slope, aspect)
    cell_height             : cy.float    = spread_inputs.cell_height
    cell_width              : cy.float    = spread_inputs.cell_width
    cell_horizontal_area    : cy.float    = cell_height * cell_width
    expected_firebrand_count: cy.float    = spot.expected_firebrand_production(fire_behavior,
                                                                               elevation_gradient,
                                                                               cell_horizontal_area,
                                                                               spot_config.firebrands_per_unit_heat)
    num_firebrands          : cy.longlong = random_generator.next_poisson(expected_firebrand_count)
    if num_firebrands > 0:
        wind_speed_10m  : cy.float               = spread_inputs.wind_speed_10m.get(t_cast, y, x)
        upwind_direction: cy.float               = spread_inputs.upwind_direction.get(t_cast, y, x)
        new_ignitions   : tuple[float, set]|None = spot.spread_firebrands(spread_inputs.fuel_model,
                                                                          spread_inputs.temperature,
                                                                          spread_inputs.fuel_moisture_dead_1hr,
                                                                          fire_type_matrix,
                                                                          (spread_inputs.rows,
                                                                           spread_inputs.cols),
                                                                          cell_height,
                                                                          cell_width,
                                                                          (t_cast, y, x),
                                                                          wind_speed_10m,
                                                                          upwind_direction,
                                                                          fire_behavior.fireline_intensity,
                                                                          fire_behavior.flame_length,
                                                                          time_of_arrival,
                                                                          random_generator,
                                                                          num_firebrands,
                                                                          spot_config)
        if new_ignitions:
            ignition_time           : float    = new_ignitions[0]
            ignited_cells           : set      = new_ignitions[1]
            concurrent_ignited_cells: set|None = spot_ignitions.get(ignition_time)
            if concurrent_ignited_cells:
                spot_ignitions[ignition_time] = set.union(ignited_cells, concurrent_ignited_cells)
            else:
                spot_ignitions[ignition_time] = ignited_cells


@cy.cfunc
@cy.exceptval(check=False)
def calc_phi_magnitude_xyz_2(phi_gradient_xy   : vec_xy,
                             elevation_gradient: vec_xy) -> cy.float:
    """
    Calculates the squared magnitude of the 3D slope-tangential phi gradient given:
    - phi_gradient_xy    :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - elevation_gradient :: (dz_dx: rise/run, dz_dy: rise/run)
    """
    (dz_dx, dz_dy)     = elevation_gradient
    (dphi_dx, dphi_dy) = phi_gradient_xy
    dphi_dz: cy.float  = vu.dot_2d(phi_gradient_xy, elevation_gradient)
    return (dphi_dx * dphi_dx + dphi_dy * dphi_dy - dphi_dz * dphi_dz / (1.0 + dz_dx * dz_dx + dz_dy * dz_dy))


# TODO: Is it faster if we save this as a top-level constant?
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def zero_partialed_wavelet() -> PartialedEllWavelet:
    return PartialedEllWavelet(
        Vh_3d = (0.0, 0.0, 0.0),
        ewc_A = 0.0,
        ewc_B = 0.0,
        ewc_C = 0.0,
    )


@cy.cfunc
@cy.exceptval(check=False)
def prepare_partialed_wavelet(heading_spread_vector: vec_xyz,
                              flanking_spread_rate : cy.float,
                              backing_spread_rate  : cy.float) -> PartialedEllWavelet:
    """
    Partially evaluates the elliptical wavelet calculation given:
    - heading_spread_vector :: (x: m/min, y: m/min, z: m/min)
    - flanking_spread_rate  :: m/min
    - backing_spread_rate   :: m/min

    The returned data is meant to be passed to function `dphi_dt_from_partialed_wavelet()`.

    NOTE: Accepting the (flanking_spread_rate, backing_spread_rate) pair is more robust than just eccentricity
          (or equivalently length_to_width_ratio) because not all fire models use elliptical wavelets
          that grow around their focus. For example, ELMFIRE defaults to something else.
    """
    heading_spread_rate: cy.float = vu.vector_magnitude_3d(heading_spread_vector)
    if heading_spread_rate > 0.0:
        heading_spread_rate_inv: cy.float = 1.0 / heading_spread_rate
        length_to_width_ratio  : cy.float = (heading_spread_rate + backing_spread_rate) / (2.0 * flanking_spread_rate)
        return PartialedEllWavelet(
            Vh_3d = heading_spread_vector,
            ewc_A = -0.5 * (heading_spread_rate - backing_spread_rate) * heading_spread_rate_inv,
            ewc_B = -(flanking_spread_rate * heading_spread_rate_inv),
            ewc_C = (length_to_width_ratio * length_to_width_ratio - 1.0),
        )
    else:
        return zero_partialed_wavelet()


@cy.cfunc
@cy.exceptval(check=False)
def wavelet_from_FireBehaviorMax(fire_behavior_max: FireBehaviorMax) -> PartialedEllWavelet:
    heading_spread_rate: cy.float = fire_behavior_max.max_spread_rate # m/min
    if heading_spread_rate > 0.0:
        heading_spread_direction: vec_xyz  = fire_behavior_max.max_spread_direction      # (x,y,z) unit vector
        heading_spread_vector   : vec_xyz  = vu.scale_3d(heading_spread_rate,
                                                         heading_spread_direction)       # (x,y,z) m/min vector
        length_to_width_ratio   : cy.float = fire_behavior_max.length_to_width_ratio     # unitless
        eccentricity            : cy.float = fire_behavior_max.eccentricity              # unitless
        backing_adjustment      : cy.float = (1.0 - eccentricity) / (1.0 + eccentricity) # unitless
        backing_spread_rate     : cy.float = heading_spread_rate * backing_adjustment    # m/min
        flanking_spread_rate    : cy.float = ((heading_spread_rate + backing_spread_rate)
                                              / (2.0 * length_to_width_ratio))           # m/min
        return prepare_partialed_wavelet(heading_spread_vector, flanking_spread_rate, backing_spread_rate)
    else:
        return zero_partialed_wavelet()


# TODO: Change local variable names
@cy.cfunc
@cy.exceptval(check=False)
def dphi_dt_from_partialed_wavelet(wavelet            : PartialedEllWavelet,
                                   phi_gradient_xy    : vec_xy,
                                   phi_magnitude_xyz_2: cy.float) -> cy.float:
    """
    Calculates the dphi/dt (a negative number in phi/min) of one elliptical wavelet given:
    - wavelet             :: PartialedEllWavelet struct
      - Vh_3d                :: (m/min, m/min, m/min)
      - ewc_A                :: unitless
      - ewc_B                :: unitless
      - ewc_C                :: unitless
    - phi_gradient_xy     :: (dphi_dx: phi/m, dphi_dy: phi/m) 2D vector on the horizontal plane
    - phi_magnitude_xyz_2 :: (phi/m)^2 squared magnitude of the 3D slope-tangential phi gradient
    """
    # Unpack vectors
    (Vx, Vy, Vz)       = wavelet.Vh_3d
    (dphi_dx, dphi_dy) = phi_gradient_xy
    # Compute intermediate values
    Vh2    : cy.float = (Vx * Vx + Vy * Vy + Vz * Vz) # NOTE: Pre-computing this doesn't appear to make it faster.
    delta  : cy.float = (Vx * dphi_dx + Vy * dphi_dy) # dot-product between Vh_3d and slope-tangential phi gradient
    # Compute dphi_dt
    return (
        wavelet.ewc_A * delta +
        wavelet.ewc_B * sqrt(
            Vh2 * phi_magnitude_xyz_2 +
            wavelet.ewc_C * (delta * delta)
        )
    )


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def phi_aware_crowning_check(phi_magnitude_xyz_2 : cy.float,
                             surface_dphi_dt     : cy.float,
                             crowning_spread_rate: cy.float) -> cy.bint:
    """
    Logically equivalent to: (surface_spread_rate > crowning_spread_rate)
    but faster to compute and robust to zero phi gradient.
    """
    return (surface_dphi_dt * surface_dphi_dt) > (crowning_spread_rate * crowning_spread_rate * phi_magnitude_xyz_2)


# NOTE: Changing this function to accept a pointer to an EllipticalInfo did not yield appreciable performance gains.
@cy.cfunc
@cy.exceptval(check=False)
def dphi_dt_from_ellipses(ellipses: EllipticalInfo, phi_gradient_xy: vec_xy) -> cy.float:
    """
    Calculates the dphi/dt (a negative number in phi/min) of the combined surface and crown elliptical wavelets.

    NOTE: The reason for computing and returning only dphi/dt is efficiency.
          Nothing else is needed in the front-propagating tight loop that iterates over tracked cells.
    """
    phi_magnitude_xyz_2: cy.float = calc_phi_magnitude_xyz_2(phi_gradient_xy, ellipses.elevation_gradient)
    surface_dphi_dt    : cy.float = dphi_dt_from_partialed_wavelet(ellipses.surface_wavelet,
                                                                   phi_gradient_xy,
                                                                   phi_magnitude_xyz_2)
    if phi_aware_crowning_check(phi_magnitude_xyz_2, surface_dphi_dt, ellipses.crowning_spread_rate):
        crown_dphi_dt: cy.float = dphi_dt_from_partialed_wavelet(ellipses.crown_wavelet,
                                                                 phi_gradient_xy,
                                                                 phi_magnitude_xyz_2)
        return min(surface_dphi_dt, crown_dphi_dt) # NOTE: dphi_dt <= 0
    else:
        return surface_dphi_dt


# TODO: Rename this constant
p_CellInputs = cy.declare(pyidx, 17) # the number of input columns


@cy.cclass
class TrackedCellsArrays:
    """
    Arrays used as on-heap supporting data structures during spread, which:
    - provide an efficient memory layout for the data involved in the Runge-Kutta passes
    - should be read-only up to the number of tracked cells
    - have a greater length in order to implement dynamic resizing

    NOTE: The *_refreshed timestamps indicate when the data was last updated for each column.

    NOTE: The motivation for cube_cache being an array of floats and not of structs
          is to enable more generic processing when reading inputs.
    """
    _array_length    : pyidx                       # power of 2, doubles each time there is a dynamic resizing
    num_tracked_cells: pyidx
    cube_cache       : cy.float[:,::1]             # shape: (_array_length, p_CellInputs)
    phi_cache        : cy.float[:,::1]             # shape: (_array_length, 9)
    sfmin_cache      : cy.pointer(FireBehaviorMin) # array of structs, caching the FireBehaviorMin for each tracked cell
    ellipse_cache    : cy.pointer(EllipticalInfo)  # array of structs (needs to be iterated over very efficiently)
    pass1_cache      : cy.pointer(Pass1CellOutput) # per-cell data produced by the 1st Runge-Kutta pass
    time_refreshed   : cy.float[17]                # an exact instant in minutes
    t_refreshed      : pyidx[17]                   # a rounded index into the space-time cubes


    def __cinit__(self,
                  time_refreshed_init: cy.float,
                  t_refreshed_init   : pyidx,
                  array_length       : pyidx = 256) -> cy.void:
        self._array_length     = array_length
        self.num_tracked_cells = 0
        self.cube_cache        = np.zeros((array_length, p_CellInputs), dtype=np.float32)
        self.phi_cache         = np.zeros((array_length, 9), dtype=np.float32)
        self.sfmin_cache       = cy.cast(cy.pointer(FireBehaviorMin),
                                         malloc(array_length * cy.sizeof(FireBehaviorMin)))
        self.ellipse_cache     = cy.cast(cy.pointer(EllipticalInfo),
                                         malloc(array_length * cy.sizeof(EllipticalInfo)))
        self.pass1_cache       = cy.cast(cy.pointer(Pass1CellOutput),
                                         malloc(array_length * cy.sizeof(Pass1CellOutput)))
        # Verify that all arrays were created
        if not (self.sfmin_cache and self.ellipse_cache and self.pass1_cache):
            # Something went wrong with malloc
            raise MemoryError()
        # Set all values in the time_refreshed and t_refreshed arrays to the passed in values
        k: pyidx
        for k in range(17):
            self.time_refreshed[k] = time_refreshed_init
            self.t_refreshed[k]    = t_refreshed_init


    @cy.cfunc
    def reset_size(self, num_tracked_cells: pyidx) -> cy.void:
        """
        Ensures that this can hold at least `num_tracked_cells`, resizing the internal arrays if necessary.
        Also updates `self.num_tracked_cells`.
        This can erase any data present in this object, so callers must make sure this information is no longer needed.
        """
        array_length: pyidx = self._array_length
        if num_tracked_cells > array_length:
            # Calculate an array_length >= num_tracked_cells
            while num_tracked_cells > array_length:
                array_length *= 2
            # Dynamically resize all internal arrays
            cube_cache   : cy.float[:,::1]             = np.zeros((array_length, p_CellInputs), dtype=np.float32)
            phi_cache    : cy.float[:,::1]             = np.zeros((array_length, 9), dtype=np.float32)
            sfmin_cache  : cy.pointer(FireBehaviorMin) = cy.cast(cy.pointer(FireBehaviorMin),
                                                                 realloc(self.sfmin_cache,
                                                                         array_length * cy.sizeof(FireBehaviorMin)))
            ellipse_cache: cy.pointer(EllipticalInfo)  = cy.cast(cy.pointer(EllipticalInfo),
                                                                 realloc(self.ellipse_cache,
                                                                         array_length * cy.sizeof(EllipticalInfo)))
            pass1_cache  : cy.pointer(Pass1CellOutput) = cy.cast(cy.pointer(Pass1CellOutput),
                                                                 realloc(self.pass1_cache,
                                                                         array_length * cy.sizeof(Pass1CellOutput)))
            # Verify that all arrays were created
            if not (sfmin_cache and ellipse_cache and pass1_cache):
                # Something went wrong with malloc
                raise MemoryError()
            else:
                self._array_length = array_length
                self.cube_cache    = cube_cache
                self.phi_cache     = phi_cache
                self.sfmin_cache   = sfmin_cache
                self.ellipse_cache = ellipse_cache
                self.pass1_cache   = pass1_cache
        # Update num_tracked_cells
        self.num_tracked_cells = num_tracked_cells


    def __dealloc__(self) -> cy.void:
        free(self.sfmin_cache)
        free(self.ellipse_cache)
        free(self.pass1_cache)


@cy.cfunc
@cy.exceptval(check=False)
def collect_phi_cache(phi_matrix: cy.float[:,::1], tca: TrackedCellsArrays) -> cy.void:
    """
    Iterates over the tracked cells and stores their phi values in a cache.
    For each tracked cell, stores a row of 9 values in `tca.phi_cache`,
    corresponding to the 'cross' of cells required to compute flux-limited gradients.
    Reads from `phi_matrix` and mutates `tca.phi_cache`.
    """
    phi_cache    : cy.float[:,::1]            = tca.phi_cache
    ellipse_cache: cy.pointer(EllipticalInfo) = tca.ellipse_cache
    i            : pyidx
    y            : pyidx
    x            : pyidx
    for i in range(tca.num_tracked_cells):
        (y, x)          = ellipse_cache[i].cell_index
        phi_cache[i, 0] = phi_matrix[2+y  , 2+x]
        phi_cache[i, 1] = phi_matrix[2+y  , 2+x-2]
        phi_cache[i, 2] = phi_matrix[2+y  , 2+x-1]
        phi_cache[i, 3] = phi_matrix[2+y  , 2+x+1]
        phi_cache[i, 4] = phi_matrix[2+y  , 2+x+2]
        phi_cache[i, 5] = phi_matrix[2+y-2, 2+x]
        phi_cache[i, 6] = phi_matrix[2+y-1, 2+x]
        phi_cache[i, 7] = phi_matrix[2+y+1, 2+x]
        phi_cache[i, 8] = phi_matrix[2+y+2, 2+x]


@cy.cfunc
@cy.exceptval(check=False)
def compare_cell_indexes(c0: coord_yx, c1: coord_yx) -> cy.int:
    """
    Lexicographic comparison of (y, x) coordinates. Returns -1, 0, or 1.
    """
    (y0, x0) = c0
    (y1, x1) = c1
    if y0 < y1:
        return -1
    elif y0 > y1:
        return 1
    # Now we know both have the same y
    elif x0 < x1:
        return -1
    elif x0 > x1:
        return 1
    else:
        return 0


# TODO: OPTIM Maybe we want to use a native array directly instead of a MemoryView.
@cy.cfunc
@cy.exceptval(check=False)
def copy_tracked_cell_data(i_old  : pyidx,
                           tca_old: TrackedCellsArrays,
                           i_new  : pyidx,
                           tca_new: TrackedCellsArrays) -> cy.void:
    # Copy cube_cache
    # NOTE: Unrolling this loop made the code 2x faster.
    tca_new.cube_cache[i_new,  0] = tca_old.cube_cache[i_old,  0]
    tca_new.cube_cache[i_new,  1] = tca_old.cube_cache[i_old,  1]
    tca_new.cube_cache[i_new,  2] = tca_old.cube_cache[i_old,  2]
    tca_new.cube_cache[i_new,  3] = tca_old.cube_cache[i_old,  3]
    tca_new.cube_cache[i_new,  4] = tca_old.cube_cache[i_old,  4]
    tca_new.cube_cache[i_new,  5] = tca_old.cube_cache[i_old,  5]
    tca_new.cube_cache[i_new,  6] = tca_old.cube_cache[i_old,  6]
    tca_new.cube_cache[i_new,  7] = tca_old.cube_cache[i_old,  7]
    tca_new.cube_cache[i_new,  8] = tca_old.cube_cache[i_old,  8]
    tca_new.cube_cache[i_new,  9] = tca_old.cube_cache[i_old,  9]
    tca_new.cube_cache[i_new, 10] = tca_old.cube_cache[i_old, 10]
    tca_new.cube_cache[i_new, 11] = tca_old.cube_cache[i_old, 11]
    tca_new.cube_cache[i_new, 12] = tca_old.cube_cache[i_old, 12]
    tca_new.cube_cache[i_new, 13] = tca_old.cube_cache[i_old, 13]
    tca_new.cube_cache[i_new, 14] = tca_old.cube_cache[i_old, 14]
    tca_new.cube_cache[i_new, 15] = tca_old.cube_cache[i_old, 15]
    tca_new.cube_cache[i_new, 16] = tca_old.cube_cache[i_old, 16]
    # Copy sfmin_cache and ellipse_cache
    # NOTE: tca_old.pass1_cache does not need to be copied over given how it will get used.
    tca_new.sfmin_cache[i_new]   = tca_old.sfmin_cache[i_old]
    tca_new.ellipse_cache[i_new] = tca_old.ellipse_cache[i_old]


# TODO: Is it faster if we save this as a top-level constant?
@cy.cfunc
@cy.inline
def inputs_name_list() -> list[str]:
    return [
        "slope",
        "aspect",
        "fuel_model",
        "canopy_cover",
        "canopy_height",
        "canopy_base_height",
        "canopy_bulk_density",
        "wind_speed_10m",
        "upwind_direction",
        "fuel_moisture_dead_1hr",
        "fuel_moisture_dead_10hr",
        "fuel_moisture_dead_100hr",
        "fuel_moisture_live_herbaceous",
        "fuel_moisture_live_woody",
        "foliar_moisture",
        "fuel_spread_adjustment",
        "weather_spread_adjustment",
    ]


@cy.cclass
class FireBehaviorSettings:
    """
    A fast-access data structure for fire behavior parameters
    to reduce the number of arguments being passed around.
    """
    max_cells_per_timestep: cy.float # CFL condition
    buffer_width          : pyidx
    use_wind_limit        : cy.bint
    surface_lw_ratio_model: str
    crown_max_lw_ratio    : cy.float
    spot_config           : dict
    cube_refresh_rates    : cy.float[17] # (min^-1) the rate at which each input column needs to be refreshed


    def __init__(self,
                 max_cells_per_timestep: float|None = 0.4,
                 buffer_width          : int|None   = 3,
                 use_wind_limit        : bool|None  = True,
                 surface_lw_ratio_model: str|None   = "behave",
                 crown_max_lw_ratio    : float|None = 1e10,
                 spot_config           : dict|None  = None,
                 cube_refresh_rates    : dict|None  = {}) -> cy.void:
        self.max_cells_per_timestep = max_cells_per_timestep
        self.buffer_width           = buffer_width
        self.use_wind_limit         = use_wind_limit
        self.surface_lw_ratio_model = surface_lw_ratio_model
        self.crown_max_lw_ratio     = crown_max_lw_ratio
        if spot_config:
            self.spot_config = SpotConfig(
                random_seed                  = spot_config.get("random_seed"),
                firebrands_per_unit_heat     = spot_config.get("firebrands_per_unit_heat", 5e-11),
                downwind_distance_mean       = spot_config.get("downwind_distance_mean", 11.7),
                fireline_intensity_exponent  = spot_config.get("fireline_intensity_exponent"),
                wind_speed_exponent          = spot_config.get("wind_speed_exponent"),
                downwind_variance_mean_ratio = spot_config.get("downwind_variance_mean_ratio"),
                crosswind_distance_stdev     = spot_config.get("crosswind_distance_stdev", 10.0),
                decay_distance               = spot_config.get("decay_distance", 500.0),
            )
        inputs_names: list = inputs_name_list()
        for k in range(17):
            self.cube_refresh_rates[k] = cube_refresh_rates[inputs_names[k]]


@cy.cfunc
@cy.exceptval(check=False)
def load_cube_cache_for_cell(spread_inputs   : SpreadInputs,
                             cell_index      : coord_yx,
                             tca             : TrackedCellsArrays,
                             i               : pyidx) -> cy.void: # NOTE: Maybe return the CellInputs struct instead?
    """
    Reads variables from input ISpaceTimeCubes and saves them by mutating `tca.cube_cache`.
    """
    (y, x)                      = cell_index
    tr        : pyidx[17]       = tca.t_refreshed
    cube_cache: cy.float[:,::1] = tca.cube_cache

    # Topography, Fuel Model, and Vegetation
    cube_cache[i,  0] = spread_inputs.slope.get(tr[0], y, x)               # rise/run
    cube_cache[i,  1] = spread_inputs.aspect.get(tr[1], y, x)              # degrees clockwise from North
    cube_cache[i,  2] = spread_inputs.fuel_model.get(tr[2], y, x)          # integer index in fm.fuel_model_table
    cube_cache[i,  3] = spread_inputs.canopy_cover.get(tr[3], y, x)        # 0-1
    cube_cache[i,  4] = spread_inputs.canopy_height.get(tr[4], y, x)       # m
    cube_cache[i,  5] = spread_inputs.canopy_base_height.get(tr[5], y, x)  # m
    cube_cache[i,  6] = spread_inputs.canopy_bulk_density.get(tr[6], y, x) # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    cube_cache[i,  7] = spread_inputs.wind_speed_10m.get(tr[7], y, x)                 # km/hr
    cube_cache[i,  8] = spread_inputs.upwind_direction.get(tr[8], y, x)               # degrees clockwise from North
    cube_cache[i,  9] = spread_inputs.fuel_moisture_dead_1hr.get(tr[9], y, x)         # kg moisture/kg ovendry weight
    cube_cache[i, 10] = spread_inputs.fuel_moisture_dead_10hr.get(tr[10], y, x)       # kg moisture/kg ovendry weight
    cube_cache[i, 11] = spread_inputs.fuel_moisture_dead_100hr.get(tr[11], y, x)      # kg moisture/kg ovendry weight
    cube_cache[i, 12] = spread_inputs.fuel_moisture_live_herbaceous.get(tr[12], y, x) # kg moisture/kg ovendry weight
    cube_cache[i, 13] = spread_inputs.fuel_moisture_live_woody.get(tr[13], y, x)      # kg moisture/kg ovendry weight
    cube_cache[i, 14] = spread_inputs.foliar_moisture.get(tr[14], y, x)               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    cube_cache[i, 15] = (spread_inputs.fuel_spread_adjustment.get(tr[15], y, x)
                         if spread_inputs.fuel_spread_adjustment is not None
                         else 1.0)                                                       # float >= 0.0
    cube_cache[i, 16] = (spread_inputs.weather_spread_adjustment.get(tr[16], y, x)
                         if spread_inputs.weather_spread_adjustment is not None
                         else 1.0)                                                       # float >= 0.0


@cy.cfunc
@cy.inline
def list_float_input_cubes(spread_inputs: SpreadInputs) -> list[ISpaceTimeCube]:
    return [
        spread_inputs.slope,
        spread_inputs.aspect,
        spread_inputs.fuel_model,
        spread_inputs.canopy_cover,
        spread_inputs.canopy_height,
        spread_inputs.canopy_base_height,
        spread_inputs.canopy_bulk_density,
        spread_inputs.wind_speed_10m,
        spread_inputs.upwind_direction,
        spread_inputs.fuel_moisture_dead_1hr,
        spread_inputs.fuel_moisture_dead_10hr,
        spread_inputs.fuel_moisture_dead_100hr,
        spread_inputs.fuel_moisture_live_herbaceous,
        spread_inputs.fuel_moisture_live_woody,
        spread_inputs.foliar_moisture,
        spread_inputs.fuel_spread_adjustment,
        spread_inputs.weather_spread_adjustment,
    ]


@cy.cfunc
def default_cube_refresh_rates(band_duration: cy.float) -> dict:
    refresh_rate: cy.float = 1.0 / band_duration
    return {
        # Non-weather inputs default to a refresh rate of 0.0 (never refreshed).
        "slope"              : 0.0,
        "aspect"             : 0.0,
        "fuel_model"         : 0.0,
        "canopy_cover"       : 0.0,
        "canopy_height"      : 0.0,
        "canopy_base_height" : 0.0,
        "canopy_bulk_density": 0.0,
        # Weather inputs default to have the same refresh rate as the base resolution of inputs.
        "wind_speed_10m"               : refresh_rate,
        "upwind_direction"             : refresh_rate,
        "fuel_moisture_dead_1hr"       : refresh_rate,
        "fuel_moisture_dead_10hr"      : refresh_rate,
        "fuel_moisture_dead_100hr"     : refresh_rate,
        "fuel_moisture_live_herbaceous": refresh_rate,
        "fuel_moisture_live_woody"     : refresh_rate,
        "foliar_moisture"              : refresh_rate,
        "fuel_spread_adjustment"       : refresh_rate,
        "weather_spread_adjustment"    : refresh_rate,
    }


recompute_levels_list: list = [
    100, # slope
    100, # aspect
    100, # fuel_model
    100, # canopy_cover
    100, # canopy_height
    100, # canopy_base_height
    100, # canopy_bulk_density
     10, # wind_speed_10m
     10, # upwind_direction
    100, # fuel_moisture_dead_1hr
    100, # fuel_moisture_dead_10hr
    100, # fuel_moisture_dead_100hr
    100, # fuel_moisture_live_herbaceous
    100, # fuel_moisture_live_woody
    100, # foliar_moisture
    100, # fuel_spread_adjustment
    100  # weather_spread_adjustment
]


# TODO: Make this more efficient by replacing the list with an array of integers.
@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def recompute_level_for_input(input_k: pyidx) -> cy.uint:
    return recompute_levels_list[input_k]


@cy.cfunc
@cy.exceptval(check=False)
def refresh_inputs_if_needed(spread_inputs: SpreadInputs,
                             fb_opts      : FireBehaviorSettings,
                             tca          : TrackedCellsArrays,
                             present_time : cy.float) -> cy.uint:
    """
    Refreshes the data input columns and refresh timestamps if needed.
    Mutates `tca` and returns an integer indicating which downstream computations need to be recomputed.
    """
    stc_list       : list[ISpaceTimeCube]|None = None
    recompute_level: cy.uint                   = 0
    k              : pyidx
    for k in range(p_CellInputs):
        refresh_rate: cy.float = fb_opts.cube_refresh_rates[k]
        # Check whether ISpaceTimeCube k needs to be refreshed
        if refresh_rate * (present_time - tca.time_refreshed[k]) > 1.0:
            # Retrieve the stc_list only the first time a ISpaceTimeCube needs to be refreshed
            stc_list = list_float_input_cubes(spread_inputs) if (stc_list is None) else stc_list
            # Extract ISpaceTimeCube k
            space_time_cube: ISpaceTimeCube = stc_list[k]
            # Calculate and store the new time_refreshed and t_refreshed values
            #
            # NOTE: The refresh time is always an integer multiple of
            #       the refresh interval. We might want to change this.
            refresh_interval  : cy.float = 1.0 / refresh_rate
            time_refreshed_new: cy.float = (present_time // refresh_interval) * refresh_interval
            t_refreshed_new   : pyidx    = int(floor(time_refreshed_new / spread_inputs.band_duration))
            tca.time_refreshed[k]        = time_refreshed_new
            tca.t_refreshed[k]           = t_refreshed_new
            # Update the cube_cache array with the latest ISpaceTimeCube values for all tracked cells
            cube_cache: cy.float[:,::1] = tca.cube_cache
            i         : pyidx
            y         : pyidx
            x         : pyidx
            for i in range(tca.num_tracked_cells):
                (y, x)           = tca.ellipse_cache[i].cell_index
                cube_cache[i, k] = space_time_cube.get(t_refreshed_new, y, x)
            # Reset the recompute_level
            recompute_level = max(recompute_level, recompute_level_for_input(k))
    return recompute_level


@cy.cfunc
@cy.exceptval(check=False)
def load_saved_CellInputs(cube_cache: cy.float[:,::1], i: pyidx) -> CellInputs:
    """
    Loads the CellInputs struct by reading the data saved in the the cube_cache array.
    """
    # Topography, Fuel Model, and Vegetation
    slope              : cy.float = cube_cache[i, 0] # rise/run
    aspect             : cy.float = cube_cache[i, 1] # degrees clockwise from North
    fuel_model_number  : cy.float = cube_cache[i, 2] # integer index in fm.fuel_model_table
    canopy_cover       : cy.float = cube_cache[i, 3] # 0-1
    canopy_height      : cy.float = cube_cache[i, 4] # m
    canopy_base_height : cy.float = cube_cache[i, 5] # m
    canopy_bulk_density: cy.float = cube_cache[i, 6] # kg/m^3

    # Wind, Surface Moisture, and Foliar Moisture
    wind_speed_10m               : cy.float = cube_cache[i, 7]  # km/hr
    upwind_direction             : cy.float = cube_cache[i, 8]  # degrees clockwise from North
    fuel_moisture_dead_1hr       : cy.float = cube_cache[i, 9]  # kg moisture/kg ovendry weight
    fuel_moisture_dead_10hr      : cy.float = cube_cache[i, 10] # kg moisture/kg ovendry weight
    fuel_moisture_dead_100hr     : cy.float = cube_cache[i, 11] # kg moisture/kg ovendry weight
    fuel_moisture_live_herbaceous: cy.float = cube_cache[i, 12] # kg moisture/kg ovendry weight
    fuel_moisture_live_woody     : cy.float = cube_cache[i, 13] # kg moisture/kg ovendry weight
    foliar_moisture              : cy.float = cube_cache[i, 14] # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    fuel_spread_adjustment   : cy.float = cube_cache[i, 15] # float >= 0.0
    weather_spread_adjustment: cy.float = cube_cache[i, 16] # float >= 0.0

    return CellInputs(
        slope                         = slope,
        aspect                        = aspect,
        fuel_model_number             = fuel_model_number,
        canopy_cover                  = canopy_cover,
        canopy_height                 = canopy_height,
        canopy_base_height            = canopy_base_height,
        canopy_bulk_density           = canopy_bulk_density,
        wind_speed_10m                = wind_speed_10m,
        upwind_direction              = upwind_direction,
        fuel_moisture_dead_1hr        = fuel_moisture_dead_1hr,
        fuel_moisture_dead_10hr       = fuel_moisture_dead_10hr,
        fuel_moisture_dead_100hr      = fuel_moisture_dead_100hr,
        fuel_moisture_live_herbaceous = fuel_moisture_live_herbaceous,
        fuel_moisture_live_woody      = fuel_moisture_live_woody,
        foliar_moisture               = foliar_moisture,
        fuel_spread_adjustment        = fuel_spread_adjustment,
        weather_spread_adjustment     = weather_spread_adjustment,
    )


@cy.cfunc
@cy.exceptval(check=False)
def resolve_surface_no_wind_no_slope_behavior(cell_inputs: CellInputs, fuel_model: FuelModel) -> FireBehaviorMin:
    """
    Computes the no-wind/no-slope surface fire behavior for a single cell.
    """
    # Pack surface fuel moisture values into a tuple
    M_f: fclaarr = (cell_inputs.fuel_moisture_dead_1hr,
                    cell_inputs.fuel_moisture_dead_10hr,
                    cell_inputs.fuel_moisture_dead_100hr,
                    0.0, # fuel_moisture_dead_herbaceous
                    cell_inputs.fuel_moisture_live_herbaceous,
                    cell_inputs.fuel_moisture_live_woody)

    # Apply fuel moisture to fuel model
    moisturized_fuel_model: FuelModel = fm.moisturize(fuel_model, M_f)

    # Combine the fuel and weather spread rate adjustments
    spread_rate_adjustment: cy.float = cell_inputs.fuel_spread_adjustment * cell_inputs.weather_spread_adjustment

    # Calculate the no-wind/no-slope surface fire behavior
    return sf.calc_surface_fire_behavior_no_wind_no_slope(moisturized_fuel_model, spread_rate_adjustment)


@cy.cfunc
@cy.exceptval(check=False)
def resolve_surface_max_behavior(fb_opts         : FireBehaviorSettings,
                                 cell_inputs     : CellInputs,
                                 fuel_model      : FuelModel,
                                 surface_fire_min: FireBehaviorMin) -> FireBehaviorMax:
    # Convert from 10m wind speed to 20ft wind speed
    wind_speed_20ft: cy.float = conv.wind_speed_10m_to_wind_speed_20ft(cell_inputs.wind_speed_10m) # km/hr

    # Convert 20ft wind speed from km/hr to m/min
    wind_speed_20ft_m_min: cy.float = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

    # Convert from 20ft wind speed to midflame wind speed in m/min
    midflame_wind_speed: cy.float = sf.calc_midflame_wind_speed(wind_speed_20ft_m_min,                   # m/min
                                                                fuel_model.delta,                        # ft
                                                                conv.m_to_ft(cell_inputs.canopy_height), # ft
                                                                cell_inputs.canopy_cover)                # 0-1

    # Calculate surface fire behavior in the direction of maximum spread
    return sf.calc_surface_fire_behavior_max(surface_fire_min,
                                             midflame_wind_speed,
                                             cell_inputs.upwind_direction,
                                             cell_inputs.slope,
                                             cell_inputs.aspect,
                                             fb_opts.use_wind_limit,
                                             fb_opts.surface_lw_ratio_model)


# TODO: OPTIM Use elevation_gradient to avoid some polar-to-cartesian conversion.
@cy.cfunc
@cy.exceptval(check=False)
def resolve_crown_max_behavior(fb_opts    : FireBehaviorSettings,
                               cell_inputs: CellInputs,
                               fuel_model : FuelModel) -> FireBehaviorMax:
    # Extract intermediate values
    heat_of_combustion          : cy.float = conv.Btu_lb_to_kJ_kg(fuel_model.h[0]) # kJ/kg
    estimated_fine_fuel_moisture: cy.float = cell_inputs.fuel_moisture_dead_1hr    # kg moisture/kg ovendry weight

    # Calculate crown fire behavior in the direction of maximum spread
    return cf.calc_crown_fire_behavior_max(cell_inputs.canopy_height,
                                           cell_inputs.canopy_base_height,
                                           cell_inputs.canopy_bulk_density,
                                           heat_of_combustion,
                                           estimated_fine_fuel_moisture,
                                           cell_inputs.wind_speed_10m,
                                           cell_inputs.upwind_direction,
                                           cell_inputs.slope,
                                           cell_inputs.aspect,
                                           fb_opts.crown_max_lw_ratio)


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def resolve_crowning_spread_rate(cell_inputs: CellInputs, surface_fire_max: FireBehaviorMax) -> cy.float:
    """
    Computes the surface spread rate at which crown fire occurs.
    """
    return cf.van_wagner_crowning_spread_rate(surface_fire_max,
                                              cell_inputs.canopy_base_height,
                                              cell_inputs.foliar_moisture)


@cy.cfunc
@cy.exceptval(check=False)
def resolve_cell_elliptical_info(fb_opts         : FireBehaviorSettings,
                                 cell_index      : coord_yx,
                                 cell_inputs     : CellInputs,
                                 fuel_model      : FuelModel,
                                 surface_fire_min: FireBehaviorMin) -> EllipticalInfo:
    # Calculate the cell's elevation gradient
    elevation_gradient: vec_xy = calc_elevation_gradient(cell_inputs.slope, cell_inputs.aspect)

    # Calculate the surface_wavelet, crown_wavelet, and crowning_spread_rate
    surface_wavelet     : PartialedEllWavelet
    crown_wavelet       : PartialedEllWavelet
    crowning_spread_rate: cy.float
    if not fuel_model.burnable:
        surface_wavelet      = zero_partialed_wavelet()
        crown_wavelet        = zero_partialed_wavelet()
        crowning_spread_rate = 1234.5 # arbitrary positive value - this threshold will never be reached.
    else:
        surface_fire_max: FireBehaviorMax = resolve_surface_max_behavior(fb_opts,
                                                                         cell_inputs,
                                                                         fuel_model,
                                                                         surface_fire_min)
        crown_fire_max  : FireBehaviorMax = resolve_crown_max_behavior(fb_opts, cell_inputs, fuel_model)
        surface_wavelet                   = wavelet_from_FireBehaviorMax(surface_fire_max)
        crown_wavelet                     = wavelet_from_FireBehaviorMax(crown_fire_max)
        crowning_spread_rate              = resolve_crowning_spread_rate(cell_inputs, surface_fire_max)

    # Build the EllipticalInfo struct
    return EllipticalInfo(
        cell_index           = cell_index,
        elevation_gradient   = elevation_gradient,
        surface_wavelet      = surface_wavelet,
        crown_wavelet        = crown_wavelet,
        crowning_spread_rate = crowning_spread_rate,
    )


@cy.cfunc
@cy.exceptval(check=False)
def refresh_caches_from_inputs_if_needed(spread_inputs: SpreadInputs,
                                         fb_opts      : FireBehaviorSettings,
                                         tca          : TrackedCellsArrays,
                                         present_time : cy.float) -> cy.void:
    """
    If required by the refresh rates, refresh inputs and recompute the necessary downstream calcs
    for each tracked cell, such as the elliptical info and the no-wind/no-slope surface fire behavior.
    Mutates `tca`.
    """
    # Update the cached cell inputs from the spread_inputs
    recompute_level  : cy.uint         = refresh_inputs_if_needed(spread_inputs, fb_opts, tca, present_time)
    cube_cache       : cy.float[:,::1] = tca.cube_cache
    cell_inputs      : CellInputs
    fuel_model_number: pyidx
    fuel_model       : FuelModel
    i                : pyidx
    if recompute_level >= 10:
        for i in range(tca.num_tracked_cells):
            # Load the cached cell inputs from the spread_inputs
            cell_inputs = load_saved_CellInputs(cube_cache, i)
            # Load the fuel model
            fuel_model_number = cy.cast(pyidx, cell_inputs.fuel_model_number)
            fuel_model        = spread_inputs.get_fm_struct(fuel_model_number)
            # Recalculate the no-wind/no-slope surface fire behavior for all tracked cells if recompute_level >= 100
            if recompute_level >= 100:
                tca.sfmin_cache[i] = resolve_surface_no_wind_no_slope_behavior(cell_inputs, fuel_model)
            # Recalculate the surface_wavelet, crown_wavelet, and crowning_spread_rate for all tracked cells
            tca.ellipse_cache[i] = resolve_cell_elliptical_info(fb_opts,
                                                                tca.ellipse_cache[i].cell_index,
                                                                cell_inputs,
                                                                fuel_model,
                                                                tca.sfmin_cache[i])


@cy.cfunc
@cy.exceptval(check=False)
def resolve_combined_spread_behavior(spread_inputs        : SpreadInputs,
                                     fb_opts              : FireBehaviorSettings,
                                     space_time_coordinate: coord_tyx,
                                     phi_gradient_xy      : vec_xy) -> SpreadBehavior:
    """
    Similar to resolve_cell_elliptical_info, but does a more exhaustive computation
    and returns a SpreadBehavior struct.
    """
    # Project phi_gradient_xy onto the slope-tangential plane as a 3D (x,y,z) vector
    cell_inputs       : CellInputs = lookup_cell_inputs(spread_inputs, space_time_coordinate)
    elevation_gradient: vec_xy     = calc_elevation_gradient(cell_inputs.slope, cell_inputs.aspect)
    phi_gradient_xyz  : vec_xyz    = calc_phi_gradient_on_slope(phi_gradient_xy, elevation_gradient)
    phi_magnitude     : cy.float   = vu.vector_magnitude_3d(phi_gradient_xyz)
    # Load the fuel model
    fuel_model_number: pyidx     = cy.cast(pyidx, cell_inputs.fuel_model_number)
    fuel_model       : FuelModel = spread_inputs.get_fm_struct(fuel_model_number)
    if phi_magnitude == 0.0 or not fuel_model.burnable:
        # This location is not on the fire perimeter and/or is not burnable
        return unburned_SpreadBehavior(elevation_gradient, phi_gradient_xyz)
    else:
        # This location is on the fire perimeter and is burnable
        surface_fire_min   : FireBehaviorMin = resolve_surface_no_wind_no_slope_behavior(cell_inputs, fuel_model)
        surface_fire_max   : FireBehaviorMax = resolve_surface_max_behavior(fb_opts,
                                                                            cell_inputs,
                                                                            fuel_model,
                                                                            surface_fire_min)
        surface_fire_normal: SpreadBehavior  = calc_fireline_normal_behavior(surface_fire_max, phi_gradient_xyz)
        # Check whether a crown fire occurs
        crowning_spread_rate: cy.float = resolve_crowning_spread_rate(cell_inputs, surface_fire_max)
        if (surface_fire_normal.spread_rate <= crowning_spread_rate):
            return surface_fire_normal
        else:
            crown_fire_max      : FireBehaviorMax = resolve_crown_max_behavior(fb_opts, cell_inputs, fuel_model)
            crown_fire_normal   : SpreadBehavior  = calc_fireline_normal_behavior(crown_fire_max, phi_gradient_xyz)
            combined_fire_normal: SpreadBehavior  = cf.calc_combined_fire_behavior(surface_fire_normal,
                                                                                   crown_fire_normal)
            return combined_fire_normal


@cy.cfunc
@cy.exceptval(check=False)
def load_tracked_cell_data(spread_inputs: SpreadInputs,
                           fb_opts      : FireBehaviorSettings,
                           cell_index   : coord_yx,
                           tca          : TrackedCellsArrays,
                           i            : pyidx) -> cy.void:
    # Read data from spread_inputs and save it in tca.cube_cache
    load_cube_cache_for_cell(spread_inputs, cell_index, tca, i)
    # Load the cached cell inputs from the spread_inputs
    cell_inputs: CellInputs = load_saved_CellInputs(tca.cube_cache, i)
    # Load the fuel model
    fuel_model_number: pyidx     = cy.cast(pyidx, cell_inputs.fuel_model_number)
    fuel_model       : FuelModel = spread_inputs.get_fm_struct(fuel_model_number)
    # Recalculate the no-wind/no-slope surface fire behavior
    surface_fire_min: FireBehaviorMin = resolve_surface_no_wind_no_slope_behavior(cell_inputs, fuel_model)
    tca.sfmin_cache[i]                = surface_fire_min
    # Recalculate the surface_wavelet, crown_wavelet, and crowning_spread_rate
    tca.ellipse_cache[i] = resolve_cell_elliptical_info(fb_opts, cell_index, cell_inputs, fuel_model, surface_fire_min)


@cy.cfunc
@cy.exceptval(check=False)
def sync_tracked_cells_arrays(spread_inputs: SpreadInputs,
                              fb_opts      : FireBehaviorSettings,
                              tracked_cells: nbt.NarrowBandTracker,
                              tca_old      : TrackedCellsArrays,
                              tca_new      : TrackedCellsArrays) -> cy.void:
    """
    Mutates `tca_new` so that it covers the same set of cells as `tracked_cells`,
    copying data from `tca_old` where possible, and otherwise loading new data from `spread_inputs`.
    """
    # Copy time_refreshed and t_refreshed arrays from tca_old to tca_new
    tca_new.reset_size(tracked_cells.num_tracked_cells)
    tca_new.time_refreshed = tca_old.time_refreshed
    tca_new.t_refreshed    = tca_old.t_refreshed
    # Set up loop variables
    cell_old     : coord_yx = (0, 0)
    cell_new     : coord_yx
    i_old        : pyidx    = 0
    i_new        : pyidx    = 0
    exhausted_old: cy.bint  = i_old >= tca_old.num_tracked_cells
    if not(exhausted_old):
        cell_old = tca_old.ellipse_cache[i_old].cell_index
    # NOTE: This loop uses the fact that tca_old is sorted consistently with compare_cell_indexes().
    ys_list: list = tracked_cells.ys_list
    if ys_list is not None:
        s: SortedDict
        for s in ys_list:
            if s is not None:
                segment: nbt.CellsCountSegment
                for segment in s.values():
                    k: pyidx
                    y: pyidx = segment.y
                    segment_counts: cy.ushort[16] = segment.counts
                    for k in range(16):
                        if (segment_counts[k] > 0):
                            # NOTE: The `for` and `if` code above is essentially just looping over the tracked cells.
                            #       This is ugly but faster than using an Iterator pattern.
                            x: pyidx = segment.x0 + k
                            cell_new = (y, x)
                            while not(exhausted_old) and compare_cell_indexes(cell_old, cell_new) < 0:
                                # cell_old is no longer tracked: just move forward.
                                i_old        += 1
                                exhausted_old = i_old >= tca_old.num_tracked_cells
                                if not(exhausted_old):
                                    cell_old = tca_old.ellipse_cache[i_old].cell_index
                            if not(exhausted_old) and (compare_cell_indexes(cell_old, cell_new) == 0):
                                # cell_new was already tracked: copy the data.
                                copy_tracked_cell_data(i_old, tca_old, i_new, tca_new)
                            else:
                                # cell_new was not in tca_old
                                load_tracked_cell_data(spread_inputs, fb_opts, cell_new, tca_new, i_new)
                            i_new += 1


@cy.cfunc
@cy.exceptval(check=False)
def runge_kutta_pass1(max_cells_per_timestep: cy.float,
                      dy                    : cy.float,
                      dx                    : cy.float,
                      max_timestep          : cy.float,
                      tca                   : TrackedCellsArrays) -> cy.float:
    """
    1st Runge-Kutta loop over elliptical dimensions, which:
    1. Resolves dt from the CFL condition
    2. Saves a Pass1CellOutput struct for each cell

    Returns the resolved `dt` and mutates `tca.pass1_cache`.
    Reads only `tca.ellipse_cache`.
    """
    ellipse_cache: cy.pointer[EllipticalInfo]  = tca.ellipse_cache
    pass1_cache  : cy.pointer[Pass1CellOutput] = tca.pass1_cache
    # The following will be useful to compute dt based on the CFL constraint.
    # It is more convenient to first compute dt_inv, the reciprocal of dt, and then
    # dt_inv = 0 represents an infinite dt. We will later enforce that dt <= max_timestep.
    dt_inv: cy.float = 0.0
    C_dx  : cy.float = max_cells_per_timestep * dx
    C_dy  : cy.float = max_cells_per_timestep * dy
    # Now looping over tracked cells:
    phi_cache: cy.float[:,::1] = tca.phi_cache
    dx_inv   : cy.float        = 1.0 / dx
    dy_inv   : cy.float        = 1.0 / dy
    for i in range(tca.num_tracked_cells):
        ellipses  : EllipticalInfo = ellipse_cache[i]
        cell_index: coord_yx       = ellipses.cell_index
        # Calculate the non-flux-limited spatial phi gradient
        dphi_dx           : cy.float = (phi_cache[i, 3] - phi_cache[i, 2]) * dx_inv / 2.0
        dphi_dy           : cy.float = (phi_cache[i, 7] - phi_cache[i, 6]) * dy_inv / 2.0
        phi_gradient_xy   : vec_xy   = (dphi_dx, dphi_dy)
        dphi_norm2        : cy.float = (dphi_dx * dphi_dx) + (dphi_dy * dphi_dy)
        dphi_dt_flim_pass1: cy.float
        if dphi_norm2 > 0.0: # Most common case
            # Calculate the flux-limited spatial phi gradient
            dphi_dx_flim: cy.float = dx_inv * calc_dphi_flim_x(phi_cache[i, 0],
                                                               phi_cache[i, 1],
                                                               phi_cache[i, 2],
                                                               phi_cache[i, 3],
                                                               phi_cache[i, 4])
            dphi_dy_flim: cy.float = dy_inv * calc_dphi_flim_y(phi_cache[i, 0],
                                                               phi_cache[i, 5],
                                                               phi_cache[i, 6],
                                                               phi_cache[i, 7],
                                                               phi_cache[i, 8])
            # Calculate the flux-limited temporal phi gradient
            dphi_dt           : cy.float = dphi_dt_from_ellipses(ellipses, phi_gradient_xy)
            dphi_dt_correction: cy.float = (dphi_dx * dphi_dx_flim + dphi_dy * dphi_dy_flim) / dphi_norm2
            dphi_dt_flim_pass1           = (dphi_dt * dphi_dt_correction)
            # Checking the CFL condition and updating dt_inv if needed (which will be rare).
            # The code is written in this way to be fast, but it's not trivial that it's correct; proof below.
            # The CFL constraint is defined as the requirement that |Ux*dt| <= C*dx and |Uy*dt| <= C*dy,
            # in which U := (Ux, Uy) is the front-normal spread rate vector in the horizontal plane,
            # and C := max_cells_per_timestep.
            # Recall that we could express U as follows: U: vec_xy = scale_2d(-dphi_dt/dphi_norm2, dphi),
            # which follows from the facts that dphi_dt = -dot2d(U, dphi) and that U is by definition
            # positively proportional to dphi.
            # In particular, Ux = -dphi_dx * dphi_dt / dphi_norm2.
            # Our constraint (from Ux) thus becomes:
            # |dt * dphi_dx * dphi_dt / dphi_norm2| <= C * dx
            # Recalling that dt_inv := 1/dt and rearranging yields:
            # dt_inv * (C * dx) * dphi_norm2 >= |dphi_dx * dphi_dt|
            if (dt_inv * (dphi_norm2 * C_dx) < abs(dphi_dt * dphi_dx)): # dt is too large given Ux
                dt_inv = abs(dphi_dt * dphi_dx) / (dphi_norm2 * C_dx)
            # And similarly for Uy:
            if (dt_inv * (dphi_norm2 * C_dy) < abs(dphi_dt * dphi_dy)): # dt is too large given Uy
                dt_inv = abs(dphi_dt * dphi_dy) / (dphi_norm2 * C_dy)
        else:
            dphi_dt_flim_pass1 = 0.0
        # Store the old and new phi values for each cell in pass1_cache
        pass1_cache[i] = Pass1CellOutput(
            cell_index      = cell_index,
            phi_gradient_xy = phi_gradient_xy,
            dphi_dt_flim    = dphi_dt_flim_pass1,
            phi_old         = phi_cache[i, 0],
        )
    # Calculate the CFL-constrained timestep dt
    dt_inv       = max(dt_inv, 1.0 / max_timestep) # (dt <= max_timestep) iff (dt_inv >= 1/max_timestep)
    dt: cy.float = 1.0 / dt_inv
    return dt


@cy.cfunc
@cy.exceptval(check=False)
def update_phi_star(tca            : TrackedCellsArrays,
                    dt             : cy.float,
                    phi_star_matrix: cy.float[:,::1]) -> cy.void:
    """
    Mutates phi_star_matrix, by using the dt and dphi/dt computed in the 1st Runge-Kutta pass.
    To be called between Runge-Kutta passes.
    """
    pass1_cache: cy.pointer[Pass1CellOutput] = tca.pass1_cache
    i          : pyidx
    for i in range(tca.num_tracked_cells):
        pass1output: Pass1CellOutput = pass1_cache[i]
        cell_index : coord_yx        = pass1output.cell_index
        y          : pyidx           = cell_index[0]
        x          : pyidx           = cell_index[1]
        phi_star_matrix[2+y, 2+x]    = pass1output.phi_old + (dt * pass1output.dphi_dt_flim)


# NOTE: Using an Extension Type here instead of a struct because it's
#       convenient to store in Python data structures like lists and dicts.
@cy.cclass
class BurnedCellInfo:
    """
    This data structure simply records information about a burned cell.
    """
    cell_index     : coord_yx
    time_of_arrival: cy.float
    phi_gradient_xy: vec_xy
    from_spotting  : cy.bint # whether spotting is what caused the cell to ignite


@cy.cfunc
def new_BurnedCellInfo(cell_index     : coord_yx,
                       time_of_arrival: cy.float,
                       phi_gradient_xy: vec_xy,
                       from_spotting  : cy.bint) -> BurnedCellInfo:
    ret: BurnedCellInfo = BurnedCellInfo()
    ret.cell_index      = cell_index
    ret.time_of_arrival = time_of_arrival
    ret.phi_gradient_xy = phi_gradient_xy
    ret.from_spotting   = from_spotting
    return ret


@cy.cfunc
def runge_kutta_pass2(dy                : cy.float,
                      dx                : cy.float,
                      start_time        : cy.float,
                      dt                : cy.float,
                      tca               : TrackedCellsArrays,
                      phi_matrix        : cy.float[:,::1]) -> list[BurnedCellInfo]:
    """
    2nd Runge-Kutta loop, which:
    1. Updates phi_matrix
    2. Identifies cells that have just burned and returns them in a list
    Reads from `tca` and `phi_star_cache`, and mutates `phi_matrix`.
    """
    dx_inv        : cy.float                    = 1.0 / dx
    dy_inv        : cy.float                    = 1.0 / dy
    ellipse_cache : cy.pointer[EllipticalInfo]  = tca.ellipse_cache
    pass1_cache   : cy.pointer[Pass1CellOutput] = tca.pass1_cache
    phi_star_cache: cy.float[:,::1]             = tca.phi_cache
    burned_cells  : list[BurnedCellInfo]        = []
    i             : pyidx
    for i in range(tca.num_tracked_cells):
        ellipses  : EllipticalInfo = ellipse_cache[i]
        cell_index: coord_yx       = ellipses.cell_index
        y         : pyidx          = cell_index[0]
        x         : pyidx          = cell_index[1]
        # Calculate the non-flux-limited spatial phi gradient
        dphi_dx           : cy.float = (phi_star_cache[i, 3] - phi_star_cache[i, 2]) * dx_inv / 2.0
        dphi_dy           : cy.float = (phi_star_cache[i, 7] - phi_star_cache[i, 6]) * dy_inv / 2.0
        phi_gradient_xy   : vec_xy   = (dphi_dx, dphi_dy)
        dphi_norm2        : cy.float = (dphi_dx * dphi_dx) + (dphi_dy * dphi_dy)
        dphi_dt_flim_pass2: cy.float
        if dphi_norm2 > 0.0: # Most common case
            # Calculate the flux-limited spatial phi gradient
            dphi_dx_flim: cy.float = dx_inv * calc_dphi_flim_x(phi_star_cache[i, 0],
                                                               phi_star_cache[i, 1],
                                                               phi_star_cache[i, 2],
                                                               phi_star_cache[i, 3],
                                                               phi_star_cache[i, 4])
            dphi_dy_flim: cy.float = dy_inv * calc_dphi_flim_y(phi_star_cache[i, 0],
                                                               phi_star_cache[i, 5],
                                                               phi_star_cache[i, 6],
                                                               phi_star_cache[i, 7],
                                                               phi_star_cache[i, 8])
            # Calculate the flux-limited temporal phi gradient
            dphi_dt           : cy.float = dphi_dt_from_ellipses(ellipses, phi_gradient_xy)
            dphi_dt_correction: cy.float = (dphi_dx * dphi_dx_flim + dphi_dy * dphi_dy_flim) / dphi_norm2
            dphi_dt_flim_pass2           = (dphi_dt * dphi_dt_correction)
        else:
            dphi_dt_flim_pass2 = 0.0
        # Combine the flux-limited temporal phi gradients from both Runge-Kutta passes
        dphi_dt_flim_pass1: cy.float = pass1_cache[i].dphi_dt_flim
        phi_old           : cy.float = pass1_cache[i].phi_old
        phi_new           : cy.float = phi_old + 0.5 * (dphi_dt_flim_pass1 + dphi_dt_flim_pass2) * dt
        phi_matrix[2+y, 2+x]         = phi_new
        # Check whether this cell has just burned, and add it to the burned_cells list if so.
        # NOTE: Phi can only ever decrease, and cells with negative phi are on fire.
        #       Therefore, if phi_old and phi_new are of opposite signs, the cell has just burned.
        if (phi_old * phi_new) < 0.0:
            time_of_arrival: cy.float = start_time + dt * phi_old / (phi_old - phi_new)
            # FIXME: Here we set phi_gradient_xy to be the phi gradient from the 1st pass.
            #        However, to be consistent with the time_of_arrival, we might want to
            #        average this with the phi_gradient_xy from the 2nd pass.
            burned_cells.append(
                new_BurnedCellInfo(cell_index      = cell_index,
                                   time_of_arrival = time_of_arrival,
                                   phi_gradient_xy = pass1_cache[i].phi_gradient_xy,
                                   from_spotting   = False)
            )
    return burned_cells


@cy.cfunc
@cy.exceptval(check=False)
def process_burned_cells(spread_inputs   : SpreadInputs,
                         fb_opts         : FireBehaviorSettings,
                         spread_state    : SpreadState,
                         spot_ignitions  : SortedDict[float, set],
                         random_generator: BufferedRandGen,
                         burned_cells    : list[BurnedCellInfo]) -> cy.void:
    # Unpack spread_state
    fire_type_matrix         : cy.uchar[:,::1] = spread_state.fire_type
    spread_rate_matrix       : cy.float[:,::1] = spread_state.spread_rate
    spread_direction_matrix  : cy.float[:,::1] = spread_state.spread_direction
    fireline_intensity_matrix: cy.float[:,::1] = spread_state.fireline_intensity
    flame_length_matrix      : cy.float[:,::1] = spread_state.flame_length
    time_of_arrival_matrix   : cy.float[:,::1] = spread_state.time_of_arrival

    # Save the burned_cells fire behavior values in the spread_state matrices
    burned_cell: BurnedCellInfo
    for burned_cell in burned_cells:
        # Determine the current space_time_coordinate
        time_of_arrival      : cy.float  = burned_cell.time_of_arrival
        cell_index           : coord_yx  = burned_cell.cell_index
        t                    : pyidx     = int(time_of_arrival // spread_inputs.band_duration)
        (y, x)                           = cell_index
        space_time_coordinate: coord_tyx = (t, y, x)
        # Re-compute the spread behavior. It's OK to re-compute it because a cell burning is a relatively rare event.
        fire_behavior: SpreadBehavior = resolve_combined_spread_behavior(spread_inputs,
                                                                         fb_opts,
                                                                         space_time_coordinate,
                                                                         burned_cell.phi_gradient_xy)
        # Write to the spread_state matrices
        fire_type_matrix[y, x]          = fire_behavior.fire_type
        spread_rate_matrix[y, x]        = fire_behavior.spread_rate
        spread_direction_matrix[y, x]   = vu.spread_direction_vector_to_angle(fire_behavior.spread_direction)
        fireline_intensity_matrix[y, x] = fire_behavior.fireline_intensity
        flame_length_matrix[y, x]       = fire_behavior.flame_length
        time_of_arrival_matrix[y, x]    = time_of_arrival

        # Cast firebrands and update spot_ignitions
        if fb_opts.spot_config:
            spot_from_burned_cell(spread_inputs,
                                  fire_type_matrix,
                                  y,
                                  x,
                                  fire_behavior,
                                  time_of_arrival,
                                  random_generator,
                                  fb_opts.spot_config,
                                  spot_ignitions)


@cy.cfunc
@cy.exceptval(check=False)
def reset_phi_star(tca               : TrackedCellsArrays,
                   spot_ignited_cells: list[BurnedCellInfo],
                   phi_star_matrix   : cy.float[:,::1],
                   phi_matrix        : cy.float[:,::1]) -> cy.void:
    """
    Efficiently updates `phi_star_matrix` to match `phi_matrix`,
    by copying only the values of cells where phi has changed.
    Mutates `phi_star_matrix`, reading from `tca.pass1_cache`, `spot_ignited_cells`, and `phi_matrix`.
    """
    y: pyidx
    x: pyidx
    i: pyidx
    # First copy phi values from the tracked cells
    for i in range(tca.num_tracked_cells):
        (y, x)                    = tca.pass1_cache[i].cell_index
        phi_star_matrix[2+y, 2+x] = phi_matrix[2+y, 2+x]
    # Then copy phi values from any spot-ignited cells
    burned_cell: BurnedCellInfo
    for burned_cell in spot_ignited_cells:
        (y, x)                    = burned_cell.cell_index
        phi_star_matrix[2+y, 2+x] = phi_matrix[2+y, 2+x]


@cy.cfunc
def ignite_from_spotting(spot_ignitions : SortedDict[float, set],
                         spread_state   : SpreadState,
                         stop_time      : cy.float) -> list[BurnedCellInfo]:
    """
    Resolves the cells to be ignited by spotting in the current time step,
    returning them as a list of (y, x) tuples, and mutates `spread_state` accordingly.
    """
    ignited_cells: list[BurnedCellInfo] = []
    if len(spot_ignitions) > 0:
        phi_matrix            : cy.float[:,::1] = spread_state.phi
        time_of_arrival_matrix: cy.float[:,::1] = spread_state.time_of_arrival
        ignition_time         : cy.float
        maybe_ignited_cells   : set
        cell_index            : coord_yx
        # https://grantjenks.com/docs/sortedcontainers/sorteddict.html
        n : pyidx = spot_ignitions.bisect_left(stop_time) # number of ignition_time values smaller than stop_time
        _i: pyidx
        for _i in range(n):
            # Remove and return the smallest ignition_time
            (ignition_time, maybe_ignited_cells) = spot_ignitions.popitem(index=0)
            for cell_index in maybe_ignited_cells:
                y: pyidx = cell_index[0]
                x: pyidx = cell_index[1]
                if phi_matrix[2+y, 2+x] > 0.0: # Not burned yet
                    phi_matrix[2+y, 2+x]        = -1.0
                    time_of_arrival_matrix[y,x] = ignition_time # FIXME: REVIEW Should I use stop_time instead?
                    ignited_cells.append(new_BurnedCellInfo(
                        cell_index      = cell_index,
                        time_of_arrival = ignition_time,
                        phi_gradient_xy = (0.0, 0.0),
                        from_spotting   = True
                    ))
                    # FIXME: I need to calculate and store the fire_behavior values for these cells
    return ignited_cells


@cy.cfunc
@cy.exceptval(check=False)
def route_cell_to_diff(frontier_cells_old: set,
                       frontier_additions: set,
                       frontier_removals : set,
                       phi_matrix        : cy.float[:,::1],
                       fuel_model_cube   : ISpaceTimeCube,
                       t                 : pyidx,
                       y                 : pyidx,
                       x                 : pyidx) -> cy.void:
    """
    Determines whether the cell `(y, x)` was just added or removed from the frontier cells,
    mutating the sets `frontier_additions` and `frontier_removals` accordingly.
    Idempotent.
    """
    encoded_cell_index: object = encode_cell_index(y, x)
    if is_frontier_cell(phi_matrix, fuel_model_cube, t, y, x):
        if not (encoded_cell_index in frontier_cells_old):
            frontier_additions.add(encoded_cell_index)
    else:
        if (encoded_cell_index in frontier_cells_old):
            frontier_removals.add(encoded_cell_index)


@cy.cfunc
def diff_frontier_cells(frontier_cells_old  : set,
                        spread_ignited_cells: list[BurnedCellInfo],
                        spot_ignited_cells  : list[BurnedCellInfo],
                        phi_matrix          : cy.float[:,::1],
                        fuel_model_cube     : ISpaceTimeCube,
                        t                   : pyidx) -> tuple[set, set]:
    """
    Computes the bi-directional set difference between the old frontier cells and the new frontier cells,
    based on newly burned cells. Returns a `(cells_added, cells_dropped)` tuple of sets, containing cell indices
    encoded by `encode_cell_index`.
    """
    frontier_additions: set = set()
    frontier_removals : set = set()
    ignited_cells     : list[BurnedCellInfo]
    # NOTE: We accept two lists below instead of one to avoid paying the cost of concatenating them.
    for ignited_cells in [spread_ignited_cells, spot_ignited_cells]:
        burned_cell: BurnedCellInfo
        y          : pyidx
        x          : pyidx
        for burned_cell in ignited_cells:
            (y, x) = burned_cell.cell_index
            # NOTE: Only in the neighborhood of a burned cell can there be changes to frontier cells membership.
            # FIXME: Should we be checking the diagonal directions as well?
            route_cell_to_diff(frontier_cells_old,
                               frontier_additions,
                               frontier_removals,
                               phi_matrix,
                               fuel_model_cube,
                               t,
                               y,
                               x)
            route_cell_to_diff(frontier_cells_old,
                               frontier_additions,
                               frontier_removals,
                               phi_matrix,
                               fuel_model_cube,
                               t,
                               y-1,
                               x)
            route_cell_to_diff(frontier_cells_old,
                               frontier_additions,
                               frontier_removals,
                               phi_matrix,
                               fuel_model_cube,
                               t,
                               y+1,
                               x)
            route_cell_to_diff(frontier_cells_old,
                               frontier_additions,
                               frontier_removals,
                               phi_matrix,
                               fuel_model_cube,
                               t,
                               y,
                               x-1)
            route_cell_to_diff(frontier_cells_old,
                               frontier_additions,
                               frontier_removals,
                               phi_matrix,
                               fuel_model_cube,
                               t,
                               y,
                               x+1)
    return (frontier_additions, frontier_removals)


@cy.cfunc
def apply_frontier_diff(frontier_cells_old: set, frontier_additions: set, frontier_removals: set) -> set:
    frontier_cells_new: set = frontier_cells_old.copy()
    encoded_cell_index: object
    for encoded_cell_index in frontier_additions:
        frontier_cells_new.add(encoded_cell_index)
    for encoded_cell_index in frontier_removals:
        frontier_cells_new.discard(encoded_cell_index)
    return frontier_cells_new


@cy.cfunc
def update_tracked_cells_with_frontier_diff(tracked_cells         : nbt.NarrowBandTracker,
                                            frontier_cells_added  : set,
                                            frontier_cells_dropped: set,
                                            buffer_width          : pyidx) -> nbt.NarrowBandTracker:
    """
    TODO: Add docstring
    """
    # Increment reference counters for all cells within buffer_width of the added frontier cells
    encoded_cell_index: object
    y                 : pyidx
    x                 : pyidx
    for encoded_cell_index in frontier_cells_added:
        (y, x) = decode_cell_index(encoded_cell_index)
        nbt.inc_square_around(tracked_cells, y, x, buffer_width)
    # Decrement reference counters for all cells within buffer_width of the dropped frontier cells
    for encoded_cell_index in frontier_cells_dropped:
        (y, x) = decode_cell_index(encoded_cell_index)
        nbt.dec_square_around(tracked_cells, y, x, buffer_width)
    # Return updated tracked cells
    return tracked_cells


@cy.cfunc
def spread_one_timestep(sim_state    : dict,
                        spread_inputs: SpreadInputs,
                        fb_opts      : FireBehaviorSettings,
                        max_timestep : cy.float) -> dict:
    """
    Spreads the fire for one iteration using the eulerian level-set method, returning an updated `sim_state`.
    """
    # Unpack sim_state
    start_time      : cy.float               = sim_state["simulation_time"]
    spread_state    : SpreadState            = sim_state["spread_state"]
    phi_matrix      : cy.float[:,::1]        = spread_state.phi
    phi_star_matrix : cy.float[:,::1]        = spread_state.phi_star
    frontier_cells  : set                    = sim_state["frontier_cells"] # TODO: OPTIM Use a binary array instead?
    tracked_cells   : nbt.NarrowBandTracker  = sim_state["tracked_cells"]
    tca             : TrackedCellsArrays     = sim_state["_tracked_cells_arrays"]
    tca_old         : TrackedCellsArrays     = sim_state["_tracked_cells_arrays_old"]
    spot_ignitions  : SortedDict[float, set] = sim_state["spot_ignitions"]
    random_generator: BufferedRandGen        = sim_state["random_generator"]

    # Insert missing tracked cells
    sync_tracked_cells_arrays(spread_inputs, fb_opts, tracked_cells, tca_old, tca)
    refresh_caches_from_inputs_if_needed(spread_inputs, fb_opts, tca, start_time)
    collect_phi_cache(phi_matrix, tca)

    # Perform the first Runge-Kutta pass and save the calculated timestep dt
    dt       : cy.float = runge_kutta_pass1(fb_opts.max_cells_per_timestep,
                                            spread_inputs.cell_height,
                                            spread_inputs.cell_width,
                                            max_timestep,
                                            tca)
    stop_time: cy.float = start_time + dt

    # Now that dt is known, update phi_star_matrix
    update_phi_star(tca, dt, phi_star_matrix)
    collect_phi_cache(phi_star_matrix, tca)

    # Perform the second Runge-Kutta pass and save the newly burned cells
    burned_cells: list[BurnedCellInfo] = runge_kutta_pass2(spread_inputs.cell_height,
                                                           spread_inputs.cell_width,
                                                           start_time,
                                                           dt,
                                                           tca,
                                                           phi_matrix)

    # Process side-effects of the burned cells (outputs, etc.)
    process_burned_cells(spread_inputs, fb_opts, spread_state, spot_ignitions, random_generator, burned_cells)
    # TODO: REVIEW It is a questionable choice to call this function AFTER process_burned_cells.
    #       It may be more sensible to ignite the spotting cells first and then to process them all.
    spot_ignited_cells: list[BurnedCellInfo] = ignite_from_spotting(spot_ignitions, spread_state, stop_time)

    # Save the new phi_matrix values in phi_star_matrix
    reset_phi_star(tca, spot_ignited_cells, phi_star_matrix, phi_matrix)

    # Update the sets of frontier cells and tracked cells based on the updated phi matrix
    frontier_diff     : tuple[set, set]       = diff_frontier_cells(frontier_cells,
                                                                    burned_cells,
                                                                    spot_ignited_cells,
                                                                    phi_matrix,
                                                                    spread_inputs.fuel_model,
                                                                    int(stop_time // spread_inputs.band_duration))
    frontier_additions: set                   = frontier_diff[0]
    frontier_removals : set                   = frontier_diff[1]
    frontier_cells_new: set                   = apply_frontier_diff(frontier_cells,
                                                                    frontier_additions,
                                                                    frontier_removals)
    tracked_cells_new : nbt.NarrowBandTracker = update_tracked_cells_with_frontier_diff(tracked_cells,
                                                                                        frontier_additions,
                                                                                        frontier_removals,
                                                                                        fb_opts.buffer_width)

    # Return the updated world state
    # NOTE: We are intentionally swapping the tracked_cells_arrays
    return {
        "simulation_time"          : stop_time,
        "spread_state"             : spread_state,
        "frontier_cells"           : frontier_cells_new,
        "tracked_cells"            : tracked_cells_new,
        "_tracked_cells_arrays"    : tca_old,
        "_tracked_cells_arrays_old": tca,
        "spot_ignitions"           : spot_ignitions,
        "random_generator"         : random_generator,
    }


@cy.cfunc
def check_space_time_cubes(space_time_cubes: dict, spot_config: dict|None = None) -> cy.void:
    # Define the provided, required, and optional keys for space_time_cubes
    provided_cubes: set = set(space_time_cubes.keys())
    required_cubes: set = {
        "slope",
        "aspect",
        "fuel_model",
        "canopy_cover",
        "canopy_height",
        "canopy_base_height",
        "canopy_bulk_density",
        "wind_speed_10m",
        "upwind_direction",
        "fuel_moisture_dead_1hr",
        "fuel_moisture_dead_10hr",
        "fuel_moisture_dead_100hr",
        "fuel_moisture_live_herbaceous",
        "fuel_moisture_live_woody",
        "foliar_moisture",
    } | ({"temperature"} if spot_config else set())
    optional_cubes: set = {
        "fuel_spread_adjustment",
        "weather_spread_adjustment",
    } | ({"temperature"} if spot_config is None else set())

    # Ensure that all required_cubes are present in provided_cubes
    if not provided_cubes.issuperset(required_cubes):
        raise ValueError("The space_time_cubes dictionary is missing these required keys: "
                         + str(required_cubes.difference(provided_cubes)))

    # Ensure that only required_cubes and optional_cubes are present in provided_cubes
    if not (required_cubes | optional_cubes).issuperset(provided_cubes):
        raise ValueError("The space_time_cubes dictionary contains these unused keys: "
                         + str(provided_cubes.difference((required_cubes | optional_cubes))))

    # Ensure that all space_time_cubes values are ISpaceTimeCube objects
    cube: object
    for cube in space_time_cubes.values():
        if not(isinstance(cube, ISpaceTimeCube)):
            raise ValueError("All values in the space_time_cubes dictionary must be ISpaceTimeCube "
                             + "objects. See pyretechnics.space_time_cube for more information.")


@cy.cfunc
def check_dimensions_and_resolutions(space_time_cubes: dict,
                                     spread_state    : SpreadState,
                                     bands           : pyidx,
                                     rows            : pyidx,
                                     cols            : pyidx,
                                     band_duration   : cy.float,
                                     cell_height     : cy.float,
                                     cell_width      : cy.float) -> cy.void:
    # Ensure that all space_time_cubes have the same cube shape
    cube: ISpaceTimeCube
    for cube in space_time_cubes.values():
        if cube.shape != (bands, rows, cols):
            raise ValueError("The space_time_cubes must all share the same cube shape.")

    # Ensure that the space_time_cubes and spread_state have the same cube shape
    if spread_state.cube_shape != (bands, rows, cols):
        raise ValueError("The space_time_cubes and spread_state must share the same cube shape.")

    # Ensure that all cube resolution values are positive
    if band_duration <= 0.0 or cell_height <= 0.0 or cell_width <= 0.0:
        raise ValueError("The cube_resolution tuple may only contain positive values.")


@cy.cfunc
def check_start_and_stop_times(start_time   : cy.float,
                               max_stop_time: cy.float,
                               cube_duration: cy.float,
                               max_duration : float|None = None) -> cy.void:
    # Ensure that start_time exists within the temporal bounds of the space_time_cubes
    if not(0.0 <= start_time < cube_duration):
        raise ValueError("The start_time falls outside of the temporal bounds of the space_time_cubes.")

    # Ensure that max_duration is positive if provided
    if max_duration and max_duration <= 0.0:
        raise ValueError("The max_duration must be a positive value if provided.")

    # Ensure that the max_stop_time does not exceed the cube_duration
    if max_stop_time > cube_duration:
        raise ValueError("The start_time + max_duration exceeds the temporal bounds of the space_time_cubes.")


@cy.ccall
def spread_fire_with_phi_field(space_time_cubes      : dict[str, ISpaceTimeCube],
                               spread_state          : SpreadState,
                               cube_resolution       : tuple[cy.float, cy.float, cy.float],
                               start_time            : cy.float,
                               max_duration          : float|None             = None,
                               max_cells_per_timestep: cy.float               = 0.4,
                               buffer_width          : pyidx                  = 3,
                               use_wind_limit        : cy.bint                = True,
                               surface_lw_ratio_model: str                    = "behave",
                               crown_max_lw_ratio    : cy.float               = 1e10,
                               spot_ignitions        : dict[float, set]       = {},
                               spot_config           : dict[str, object]|None = None,
                               cube_refresh_rates    : dict[str, float]       = {}) -> dict[str, object]:
    """
    Given these inputs:
    - space_time_cubes             :: dictionary of ISpaceTimeCube objects with these cell types
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
    - spread_state                 :: SpreadState object whose spatial dimensions match the space_time_cubes
    - cube_resolution              :: tuple with these fields
      - band_duration                 :: minutes
      - cell_height                   :: meters
      - cell_width                    :: meters
    - start_time                   :: minutes (from the start of the space_time_cube's temporal origin)
    - max_duration                 :: minutes (Optional)
    - max_cells_per_timestep       :: max number of cells the fire front can travel in one timestep (Optional)
    - buffer_width                 :: Chebyshev distance from frontier cells to include in tracked cells (Optional)
    - use_wind_limit               :: boolean (Optional)
    - surface_lw_ratio_model       :: "rothermel" or "behave" (Optional)
    - crown_max_lw_ratio           :: float > 0.0 (Optional)
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
    - cube_refresh_rates           :: dictionary from input name to refresh rate in 1/min (Optional).
                                      0 means never refresh. Weather inputs default to 1/band_duration,
                                      whereas non-weather inputs default to 0.

    return a dictionary with these keys:
    - stop_time         :: minutes
    - stop_condition    :: "max duration reached" or "no burnable cells"
    - spread_state      :: SpreadState object whose spatial dimensions match the space_time_cubes
    - num_tracked_cells :: number of cells in the narrow band at stop_time
    - spot_ignitions    :: dictionary of (ignition_time -> ignited_cells) (only included when spotting is used)
    - random_generator  :: BufferedRandGen object (only included when spotting is used)
    """
    # Verify the contents of space_time_cubes
    check_space_time_cubes(space_time_cubes, spot_config)

    # Extract simulation dimensions
    fuel_model_cube: ISpaceTimeCube             = space_time_cubes["fuel_model"]
    cube_shape     : tuple[pyidx, pyidx, pyidx] = fuel_model_cube.shape
    bands          : pyidx                      = cube_shape[0]
    rows           : pyidx                      = cube_shape[1]
    cols           : pyidx                      = cube_shape[2]

    # Extract simulation resolutions
    band_duration: cy.float = cube_resolution[0]
    cell_height  : cy.float = cube_resolution[1]
    cell_width   : cy.float = cube_resolution[2]

    # Verify the simulation dimensions and resolutions
    check_dimensions_and_resolutions(space_time_cubes,
                                     spread_state,
                                     bands,
                                     rows,
                                     cols,
                                     band_duration,
                                     cell_height,
                                     cell_width)

    # Calculate the cube duration and max stop time
    cube_duration: cy.float = bands * band_duration
    max_stop_time: cy.float = start_time + max_duration if max_duration else cube_duration

    # Verify simulation start and stop times
    check_start_and_stop_times(start_time, max_stop_time, cube_duration, max_duration)

    # Identify the sets of frontier cells and tracked cells based on the phi matrix
    start_t       : pyidx                 = int(start_time // band_duration)
    frontier_cells: set                   = identify_all_frontier_cells(spread_state.phi,
                                                                        fuel_model_cube,
                                                                        start_t,
                                                                        rows,
                                                                        cols)
    tracked_cells : nbt.NarrowBandTracker = identify_tracked_cells(frontier_cells, buffer_width, rows, cols)

    # Create a BufferedRandGen object to produce random samples if spot_config is provided
    random_generator: BufferedRandGen|None = None
    if spot_config:
        random_generator = BufferedRandGen(np.random.default_rng(seed=spot_config.get("random_seed")))

    # Prepare the SpreadInputs struct
    spread_inputs: SpreadInputs = make_SpreadInputs(cube_shape, cube_resolution, space_time_cubes)

    # Prepare the FireBehaviorSettings struct
    fb_opts: FireBehaviorSettings = FireBehaviorSettings(
        max_cells_per_timestep = max_cells_per_timestep,
        buffer_width           = buffer_width,
        use_wind_limit         = use_wind_limit,
        surface_lw_ratio_model = surface_lw_ratio_model,
        crown_max_lw_ratio     = crown_max_lw_ratio,
        spot_config            = spot_config,
        cube_refresh_rates     = {**default_cube_refresh_rates(band_duration), **cube_refresh_rates},
    )

    # Prepare the sim_state dictionary
    # NOTE: We are intentionally swapping the tracked_cells_arrays.
    #       It's OK not to be in sync - spread_one_timestep will solve this.
    # TODO: Turn sim_state into a struct
    sim_state: dict = {
        "simulation_time"          : start_time,
        "spread_state"             : spread_state,
        "frontier_cells"           : frontier_cells,
        "tracked_cells"            : tracked_cells,
        "_tracked_cells_arrays"    : TrackedCellsArrays(start_time, start_t),
        "_tracked_cells_arrays_old": TrackedCellsArrays(start_time, start_t),
        "spot_ignitions"           : SortedDict(spot_ignitions), # Convert spot_ignitions into a SortedDict
        "random_generator"         : random_generator,
    }

    # FIXME: I don't think the "no burnable cells" condition can ever be met currently.
    # Spread the fire until an exit condition is reached
    remaining_time_in_simulation: cy.float = max_stop_time - start_time
    early_exit_threshold        : cy.float = 1.0 / 60.0 # 1 second
    while((remaining_time_in_simulation > early_exit_threshold)       # 1. There is still time left in the simulation
          and (nbt.nonempty_tracked_cells(tracked_cells)              # 2. There are burning cells on the grid
               or len(sim_state["spot_ignitions"]) > 0)):             # 3. There are embers waiting to catch fire
        # Spread fire one timestep
        sim_state                    = spread_one_timestep(sim_state,
                                                           spread_inputs,
                                                           fb_opts,
                                                           remaining_time_in_simulation)
        remaining_time_in_simulation = max_stop_time - sim_state["simulation_time"]

    # Determine the stop_condition
    stop_condition: str = ("max duration reached"
                           if remaining_time_in_simulation <= early_exit_threshold
                           else "no burnable cells")

    # Return the final simulation results
    return {
        "stop_time"        : sim_state["simulation_time"],
        "stop_condition"   : stop_condition,
        "spread_state"     : sim_state["spread_state"],
        "num_tracked_cells": tracked_cells.num_tracked_cells,
    } | ({
        "spot_ignitions"  : sim_state["spot_ignitions"],
        "random_generator": sim_state["random_generator"],
    } if spot_config else {})
# spread-phi-field ends here
