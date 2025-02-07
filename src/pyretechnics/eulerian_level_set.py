# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients-approx][phi-field-spatial-gradients-approx]]
# cython: profile=False
import cython
import cython as cy
import numpy as np
from sortedcontainers import SortedDict
if cython.compiled:
    from cython.cimports.cpython.mem import PyMem_Malloc, PyMem_Free # Unique to Compiled Cython
    from cython.cimports.pyretechnics.math import pi, floor, sqrt, atan
    from cython.cimports.pyretechnics.cy_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, \
        fclaarr, FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior, SpotConfig, PartialedEllWavelet
    from cython.cimports.pyretechnics.random import BufferedRandGen
    from cython.cimports.pyretechnics.space_time_cube import ISpaceTimeCube
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.vector_utils as vu
    import cython.cimports.pyretechnics.fuel_models as fm
    import cython.cimports.pyretechnics.surface_fire as sf
    import cython.cimports.pyretechnics.crown_fire as cf
    import cython.cimports.pyretechnics.spot_fire as spot
    import cython.cimports.pyretechnics.narrow_band_tracking as nbt
else:
    # TODO: Create equivalent Python functions for PyMem_Malloc, PyMem_Free
    from math import pi, floor, sqrt, atan
    from pyretechnics.py_types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx, \
        fclaarr, FuelModel, FireBehaviorMin, FireBehaviorMax, SpreadBehavior, SpotConfig, PartialedEllWavelet
    from pyretechnics.random import BufferedRandGen
    from pyretechnics.space_time_cube import ISpaceTimeCube
    import pyretechnics.conversion as conv
    import pyretechnics.vector_utils as vu
    import pyretechnics.fuel_models as fm
    import pyretechnics.surface_fire as sf
    import pyretechnics.crown_fire as cf
    import pyretechnics.spot_fire as spot
    import pyretechnics.narrow_band_tracking as nbt


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
@cy.wraparound(False)
@cy.cdivision(True)
@cy.boundscheck(False)
def calc_dphi_dx_approx(phi: cy.float[:,:], dx: cy.float, x: pyidx, y: pyidx) -> cy.float: # NOTE no longer used in tight loops.
    """
    Calculate the spatial gradient of the phi raster in the x (west->east)
    direction at grid cell (x,y) given the cell width dx.
    """
    east_x: pyidx = x + 1
    west_x: pyidx = x - 1
    return (phi[2+y][2+east_x] - phi[2+y][2+west_x]) / (2.0 * dx)


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
@cy.wraparound(False)
@cy.cdivision(True)
@cy.boundscheck(False)
def calc_dphi_dy_approx(phi: cy.float[:,:], dy: cy.float, x: pyidx, y: pyidx) -> cy.float: # NOTE no longer used in tight loops.
    """
    Calculate the spatial gradient of the phi raster in the y (south->north)
    direction at grid cell (x,y) given the cell height dy.
    """
    north_y: pyidx = y + 1
    south_y: pyidx = y - 1
    return (phi[2+north_y][2+x] - phi[2+south_y][2+x]) / (2.0 * dy)


@cy.cfunc
def calc_phi_gradient_approx(phi: cy.float[:,:], dx: cy.float, dy: cy.float, x: pyidx, y: pyidx) -> vec_xy: # NOTE no longer used in tight loops.
    """
    Calculate the spatial gradient of the phi raster at grid cell (x,y)
    given the cell width dx and the cell height dy.
    """
    dphi_dx: cy.float = calc_dphi_dx_approx(phi, dx, x, y)
    dphi_dy: cy.float = calc_dphi_dy_approx(phi, dy, x, y)
    return (dphi_dx, dphi_dy)
# phi-field-spatial-gradients-approx ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector][phi-field-normal-vector]]
# TODO: Remove unused function
@cy.cfunc
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
        return vu.as_unit_vector_2d(phi_gradient) # (n_x, n_y)
# phi-field-normal-vector ends here
# [[file:../../org/pyretechnics.org::phi-field-normal-vector-angle][phi-field-normal-vector-angle]]
# TODO: Remove unused function
@cy.cfunc
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
@cy.inline
def half_superbee_dphi_up(dphi_up: cy.float, dphi_loc: cy.float) -> cy.float:
    """
    Logically like calc_superbee_flux_limiter(), but returns a result multiplied by (0.5 * dphi_loc).
    """
    # NOTE this is more numerically stable than calc_superbee_flux_limiter().
    s_loc: cy.float = 1.0 if dphi_loc >= 0.0 else -1.0
    are_opposite_signs: cy.bint = (s_loc * dphi_up) <= 0.0
    if are_opposite_signs:
        return 0.0
    a_up: cy.float = abs(dphi_up)
    a_loc: cy.float = abs(dphi_loc)
    return s_loc * max(
        min(a_up / 2, a_loc),
        min(a_up, a_loc / 2))
# superbee-flux-limiter ends here
# [[file:../../org/pyretechnics.org::phi-field-spatial-gradients][phi-field-spatial-gradients]]
@cy.cfunc
@cy.exceptval(check=False)
def calc_dphi_flim_x(p00: cy.float, pw2: cy.float, pw1: cy.float, pe1: cy.float, pe2: cy.float) -> cy.float:
    dphi_loc: cy.float

    phi_east: cy.float
    dphi_loc = pe1 - p00
    if pe1 >= pw1:
        dphi_up: cy.float = p00 - pw1
        phi_east = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up: cy.float = pe2 - pe1
        phi_east = pe1 - half_superbee_dphi_up(dphi_up, dphi_loc)

    phi_west: cy.float
    dphi_loc = pw1 - p00
    if pe1 >= pw1:
        dphi_up: cy.float = pw2 - pw1
        phi_west = pw1 - half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up: cy.float = p00 - pe1
        phi_west = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)
    return (phi_east - phi_west)

# NOTE this is actually the same function as the previous one. But who knows, maybe we get a performance gain by differentiating code sites.
@cy.cfunc
@cy.exceptval(check=False)
def calc_dphi_flim_y(p00: cy.float, ps2: cy.float, ps1: cy.float, pn1: cy.float, pn2: cy.float) -> cy.float:
    dphi_loc: cy.float

    phi_north: cy.float
    dphi_loc = pn1 - p00
    if pn1 >= ps1:
        dphi_up: cy.float = p00 - ps1
        phi_north = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up: cy.float = pn2 - pn1
        phi_north = pn1 - half_superbee_dphi_up(dphi_up, dphi_loc)

    phi_south: cy.float
    dphi_loc = ps1 - p00
    if pn1 >= ps1:
        dphi_up: cy.float = ps2 - ps1
        phi_south = ps1 - half_superbee_dphi_up(dphi_up, dphi_loc)
    else:
        dphi_up: cy.float = p00 - pn1
        phi_south = p00 + half_superbee_dphi_up(dphi_up, dphi_loc)

    return (phi_north - phi_south)
# phi-field-spatial-gradients ends here
# [[file:../../org/pyretechnics.org::phi-east][phi-east]]
# phi-east ends here
# [[file:../../org/pyretechnics.org::phi-west][phi-west]]
# phi-west ends here
# [[file:../../org/pyretechnics.org::phi-north][phi-north]]
# phi-north ends here
# [[file:../../org/pyretechnics.org::phi-south][phi-south]]
# phi-south ends here
# [[file:../../org/pyretechnics.org::calc-fireline-normal-behavior][calc-fireline-normal-behavior]]
# TODO: Move this to pyretechnics.vector_utils and use throughout the literate program
@cy.cfunc
def calc_elevation_gradient(slope: cy.float, aspect: cy.float) -> vec_xy:
    """
    Returns the elevation gradient (dz_dx: rise/run, dz_dy: rise/run) given:
    - slope  :: rise/run
    - aspect :: degrees clockwise from North
    """
    return conv.azimuthal_to_cartesian(slope, conv.opposite_direction(aspect))


# NOTE it would be better to use a cython enum here, but that's not supported in pure python syntax...
fire_type_unburned      = cy.declare(cy.int, 0)
fire_type_surface       = cy.declare(cy.int, 1)
fire_type_crown_passive = cy.declare(cy.int, 2)
fire_type_crown_active  = cy.declare(cy.int, 3)


@cy.cfunc
def calc_phi_gradient_on_slope(phi_gradient_xy: vec_xy, elevation_gradient: vec_xy) -> vec_xyz:
    """
    Return the gradient of phi projected onto the slope-tangential plane as a 3D (x,y,z) vector (in phi/m) given:
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

    phi_magnitude = vu.vector_magnitude_3d(phi_gradient) # phi/m

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
        heading_spread_vector: vec_xyz       = vu.scale_3d(heading_spread_rate, heading_spread_direction)     # (x,y,z) m/min vector
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

        B: cy.float = vu.dot_3d(heading_spread_vector, phi_gradient)                             # phi/min
        C: cy.float = flanking_spread_rate / heading_spread_rate                              # unitless
        D: cy.float = (heading_spread_rate * phi_magnitude) ** 2.0                            # (phi/min)^2
        E: cy.float = (length_to_width_ratio ** 2.0 - 1.0) * (B ** 2.0)                       # (phi/min)^2
        dphi_dt: cy.float = -(A * B + C * sqrt(D + E))                                              # phi/min

        #============================================================================================
        # Calculate fire behavior normal to the fire perimeter
        #============================================================================================

        normal_spread_rate: cy.float        = -dphi_dt / phi_magnitude                        # m/min
        normal_direction: vec_xyz           = vu.as_unit_vector_3d(phi_gradient)                 # (x,y,z) unit vector
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
        self.spatial_resolution = (cube_resolution[2], cube_resolution[1])

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


    @cy.cfunc
    def _init_fuel_models(self):
        self.fuel_models_arr = cy.cast(
            cy.pointer(FuelModel),
            PyMem_Malloc(300 * cython.sizeof(FuelModel)))
        if not self.fuel_models_arr:
            raise MemoryError()
        f: FuelModel
        for f in fm.fuel_model_table.values():
            self.fuel_models_arr[f.number] = f

    @cy.cfunc
    def get_fm_struct(self, fm_number: pyidx) -> FuelModel:
        fuel_models_arr: cy.pointer(FuelModel) = self.fuel_models_arr
        fuel_model: FuelModel = fuel_models_arr[fm_number]
        return fuel_model


    def __dealloc__(self):
        PyMem_Free(self.fuel_models_arr)  # no-op if self.data is NULL

@cy.cfunc
def make_SpreadInputs(cube_resolution, space_time_cubes: dict) -> SpreadInputs:
    sinputs: SpreadInputs = SpreadInputs(cube_resolution,
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
    return sinputs


@cy.cfunc
def lookup_space_time_cube_float32(space_time_cube: ISpaceTimeCube, space_time_coordinate: coord_tyx) -> cy.float:
    t: pyidx = space_time_coordinate[0]
    y: pyidx = space_time_coordinate[1]
    x: pyidx = space_time_coordinate[2]
    return space_time_cube.get(t, y, x)


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
def lookup_cell_inputs(space_time_cubes: SpreadInputs, tyx: coord_tyx) -> CellInputs:
    """
    Reads the inputs for a given cell from the space-time cubes, returning a `CellInputs` struct.
    """
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

@cy.cfunc
@cy.cdivision(True)
@cy.exceptval(check=False)
def unburned_SpreadBehavior(elevation_gradient: vec_xy, dphi_st: vec_xyz) -> SpreadBehavior:
    # NOTE we're only going through these annoying calculations because we are required to return a spread_direction unit vector, which is of questionable value.
    # IMPROVEMENT We wouldn't have to go through this trouble if we simply returned a Cartesian speed vector instead, which would play more nicely with the rest of the code.
    spread_direction: vec_xyz
    if dphi_st != (0, 0, 0): # IMPROVEMENT numerical stability
        spread_direction = vu.scale_3d(1./ vu.vector_magnitude_3d(dphi_st), dphi_st)
    elif elevation_gradient != (0, 0): # IMPROVEMENT numerical stability
        slope_vector_3d: vec_xyz = vu.to_slope_plane(elevation_gradient, elevation_gradient)
        spread_direction = vu.as_unit_vector_3d(slope_vector_3d)
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
# burn-cell-toward-phi-gradient ends here
# [[file:../../org/pyretechnics.org::phi-field-perimeter-tracking][phi-field-perimeter-tracking]]
@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def opposite_phi_signs(phi_matrix: cy.float[:,:], y1: pyidx, x1: pyidx, y2: pyidx, x2: pyidx) -> cy.bint:
    """
    Return True if the phi values at cells (x1,y1) and (x2,y2) have opposite signs.
    """
    return phi_matrix[2+y1, 2+x1] * phi_matrix[2+y2, 2+x2] < 0.0


# TODO: Is it faster to build up a list or a set?
# TODO: Should we store each frontier_cells entry as a coord_yx?
@cy.cfunc
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
            north_y: pyidx = y+1
            south_y: pyidx = y-1
            east_x : pyidx = x+1
            west_x : pyidx = x-1
            if (opposite_phi_signs(phi_matrix, y, x, north_y, x) or
                opposite_phi_signs(phi_matrix, y, x, south_y, x) or
                opposite_phi_signs(phi_matrix, y, x, y, east_x) or
                opposite_phi_signs(phi_matrix, y, x, y, west_x)):
                frontier_cells.add(encode_cell_index(y, x))
    return frontier_cells

@cy.cfunc
def is_frontier_cell(phi_matrix: cy.float[:,:], y: pyidx, x: pyidx) -> cy.bint:
    north_y: pyidx = y+1
    south_y: pyidx = y-1
    east_x : pyidx = x+1
    west_x : pyidx = x-1
    return (
        opposite_phi_signs(phi_matrix, y, x, north_y, x) or
        opposite_phi_signs(phi_matrix, y, x, south_y, x) or
        opposite_phi_signs(phi_matrix, y, x, y, east_x) or
        opposite_phi_signs(phi_matrix, y, x, y, west_x)
        )


@cy.cfunc
def identify_tracked_cells(frontier_cells: set, buffer_width: pyidx, rows: pyidx, cols: pyidx) -> object:
    """
    TODO: Add docstring
    """
    tracked_cells: object = nbt.new_NarrowBandTracker(cols, rows)
    cell         : object
    y: pyidx
    x: pyidx
    for cell in frontier_cells:
        y, x = decode_cell_index(cell)
        nbt.inc_square_around(tracked_cells, y, x, buffer_width)
    return tracked_cells


# phi-field-perimeter-tracking ends here
# [[file:../../org/pyretechnics.org::spread-phi-field][spread-phi-field]]
@cy.cfunc
def spot_from_burned_cell(
        spot_config: SpotConfig,
        sinputs: SpreadInputs,
        fire_type_matrix: cy.uchar[:,:],
        random_generator: BufferedRandGen,
        y: pyidx,
        x: pyidx,
        fb: SpreadBehavior,
        toa: cy.float,
        spot_ignitions: object
    ) -> cy.void:
    """
    Schedules the future spot ignitions following from burning the given cell.
    Mutates `spot_ignitions` (and `random_generator`).
    """
    cell_horizontal_area_m2 : cy.float = sinputs.spatial_resolution[0] * sinputs.spatial_resolution[1]
    t_cast                  : pyidx    = int(toa // sinputs.band_duration)
    space_time_coordinate   : coord_tyx = (t_cast, y, x)
    slope                   : cy.float = lookup_space_time_cube_float32(sinputs.slope, space_time_coordinate)
    aspect                  : cy.float = lookup_space_time_cube_float32(sinputs.aspect, space_time_coordinate)
    elevation_gradient      : vec_xy   = calc_elevation_gradient(slope, aspect)
    firebrands_per_unit_heat: cy.float = spot_config.firebrands_per_unit_heat
    expected_firebrand_count: cy.float = spot.expected_firebrand_production(fb,
                                                                            elevation_gradient,
                                                                            cell_horizontal_area_m2,
                                                                            firebrands_per_unit_heat)
    num_firebrands: cy.long = random_generator.next_poisson(expected_firebrand_count)
    if num_firebrands > 0:
        # OPTIM we might want to hold to the SpreadInputs and look these up in there.
        wind_speed_10m: cy.float = lookup_space_time_cube_float32(sinputs.wind_speed_10m, space_time_coordinate)
        upwind_direction: cy.float = lookup_space_time_cube_float32(sinputs.upwind_direction, space_time_coordinate)
        new_ignitions: tuple[float, set]|None = spot.spread_firebrands(sinputs.fuel_model,
                                                                       sinputs.temperature,
                                                                       sinputs.fuel_moisture_dead_1hr,
                                                                       fire_type_matrix,
                                                                       (sinputs.rows, sinputs.cols),
                                                                       sinputs.spatial_resolution,
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


@cy.cfunc
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


@cy.cfunc
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
    Vh: cy.float = vu.vector_magnitude_3d(Vh_3d)
    Vh_inv: cy.float = (1.0/Vh) if Vh > 0 else 0.0
    LoW: cy.float = 0.5*(Vh + Vb)/Vf if Vf > 0 else 1.0 # Length/Width ratio
    return PartialedEllWavelet(
        # NOTE if Vh = 0, all of the following will be 0, as they should.
        Vh_3d = Vh_3d,
        ewc_A = -0.5 * (Vh - Vb) * Vh_inv, # INTRO
        ewc_B = -(Vf * Vh_inv), # INTRO
        ewc_C = (LoW*LoW - 1) # INTRO
        )


@cy.cfunc
@cy.exceptval(check=False)
def zero_partialed_wavelet() -> PartialedEllWavelet:
    Vh_3d: vec_xyz = (0, 0, 0)
    return prepare_partialed_wavelet(Vh_3d, 0, 0)


@cy.cfunc
@cy.cdivision(True)
@cy.exceptval(check=False)
def pw_from_FireBehaviorMax(fb_max: FireBehaviorMax) -> PartialedEllWavelet:
    Vh: cy.float = fb_max.max_spread_rate
    Vh_3d: vec_xyz = vu.scale_3d(Vh, fb_max.max_spread_direction)
    if Vh > 0:
        LoW: cy.float = fb_max.length_to_width_ratio
        ecc: cy.float = fb_max.eccentricity
        backing_adjustment: cy.float   = (1.0 - ecc) / (1.0 + ecc)
        Vb: cy.float  = Vh * backing_adjustment
        Vf: cy.float = (Vh + Vb) / (2.0 * LoW) # m/min
        return prepare_partialed_wavelet(Vh_3d, Vf, Vb)
    else:
        return zero_partialed_wavelet()



@cy.cfunc
@cy.exceptval(check=False)
@cy.inline
def dphi_dt_from_partialed_wavelet(
        pw: PartialedEllWavelet,
        dphi: vec_xy, # 2D horizontal gradient of phi.
        st_dphi_2: cy.float, # INTRO Squared norm (_2) of slope-tangential (st_) phi gradient (dphi)
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
    crowning_spread_rate_threshold = cy.float, # INTRO pre-computed critical threshold in surface spread rate at which crowning occurs.
    crown_wavelet = PartialedEllWavelet,
)
# NOTE a significant benefit of this architecture is that it's Rothermel-agnostic:
# EllipticalInfo could conceivably be implemented using variants of the Rothermel model.
# This can be valuable to give flexibility to users.

@cy.cfunc
@cy.exceptval(check=False)
def dphi_dt_from_elliptical(ell_i: EllipticalInfo, dphi: vec_xy) -> cy.float: # NOTE this code has to be very fast!
    """
    Calculates the dphi/dt (a negative number in phi/min) of the combined surface and crown elliptical wavelets.

    The reason for computing and returning _only_ dphi/dt is efficiency:
    nothing else is needed in the front-propagating tight loop
    that iterates over tracked cells.
    """
    # NOTE changing this function to accept a pointer to an EllipticalInfo did not yield appreciable performance gains.
    st_dphi_2: cy.float = compute_st_dphi_2(ell_i.slp_dz, dphi)
    surfc_dphi_dt: cy.float = dphi_dt_from_partialed_wavelet(ell_i.surfc_wavelet, dphi, st_dphi_2)
    csr: cy.float = ell_i.crowning_spread_rate_threshold
    csr_2: cy.float = csr * csr
    does_crown: cy.bint = (surfc_dphi_dt * surfc_dphi_dt) > csr_2 * st_dphi_2 # Logically equivalent to: (surface_spread_rate > crowning_spread_rate_threshold), but faster to compute and robust to zero phi gradient.
    if does_crown:
        crown_dphi_dt: cy.float = dphi_dt_from_partialed_wavelet(ell_i.crown_wavelet, dphi, st_dphi_2)
        return min(surfc_dphi_dt, crown_dphi_dt) # Note that dphi_dt <= 0
    else:
        return surfc_dphi_dt

Pass1CellOutput = cy.struct( # INTRO some data saved during the 1st Runge-Kutta pass.
    cell_index = coord_yx,
    dphi = vec_xy,
    dphi_dt_flim = cy.float, # Flux-limited dphi/dt (phi/min, <= 0).
    phi_old = cy.float
    )

p_CellInputs = cy.declare(pyidx, 17) # The number of input columns.

@cy.cfunc
def inputs_name_list() -> list[str]:
    return [
        'slope',
        'aspect',

        'fuel_model',

        'canopy_cover',
        'canopy_height',
        'canopy_base_height',
        'canopy_bulk_density',

        'wind_speed_10m',
        'upwind_direction',

        'fuel_moisture_dead_1hr',
        'fuel_moisture_dead_10hr',
        'fuel_moisture_dead_100hr',
        'fuel_moisture_live_herbaceous',
        'fuel_moisture_live_woody',
        'foliar_moisture',

        'fuel_spread_adjustment',
        'weather_spread_adjustment',
    ]

@cy.cclass
class TrackedCellsArrays:
    """
    Arrays used as on-heap supporting data structures during spread.
    """
    _arrays_length: pyidx # power of 2, double each time there is a dynamic resizing.
    n_tracked_cells: pyidx


    # These timestamps say when the data was last updated for each data column.
    time_refreshed: cy.float[17] # This one is an exact instant in minutes.
    t_refreshed: pyidx[17] # This one is an rounded index into the inputs datacubes.
    # REVIEW The t_ vs time_ naming is embarassingly confusion-prone, but I don't know how to help that - this confusion is already all over the place in the codebase.

    # These arrays provide an efficient memory layout for the data involved in the Runge-Kutta passes.
    # These arrays should be read only up to the number of tracked cells;
    # they have a greater length in order to implement dynamic resizing.
    float_inputs: cy.float[:, :] # Shape: (n_tracked_cells, p_CellInputs)
    # NOTE the motivation for float_inputs being an array of floats and not of structs is to enable more generic processing when reading inputs.
    sfmin_arr: cy.pointer(FireBehaviorMin) # Array of structs, caching the FireBehaviorMin for each tracked cell.
    ell_info: cy.pointer(EllipticalInfo) # Array of structs (needs to be iterated over very efficiently).
    pass1outputs: cy.pointer(Pass1CellOutput) # Per-cell data produced by the 1st Runge-Kutta pass.

    phi_values: cy.float[:,:] # Shape: (self._arrays_length, 9)

    def __cinit__(self,
            time_refreshed_init: cy.float,
            t_refreshed_init: pyidx,
            _arrays_length: pyidx = 256):
        self._arrays_length = _arrays_length
        self.n_tracked_cells = 0
        self.float_inputs = np.zeros((self._arrays_length, p_CellInputs), dtype=np.float32)
        self.sfmin_arr = cy.cast(cy.pointer(FireBehaviorMin), PyMem_Malloc(self._arrays_length * cy.sizeof(FireBehaviorMin)))
        if not self.sfmin_arr:
            raise MemoryError()
        self.ell_info = cy.cast(cy.pointer(EllipticalInfo), PyMem_Malloc(self._arrays_length * cy.sizeof(EllipticalInfo)))
        if not self.ell_info:
            raise MemoryError()
        self.pass1outputs = cy.cast(cy.pointer(Pass1CellOutput), PyMem_Malloc(self._arrays_length * cy.sizeof(Pass1CellOutput)))
        for k in range(17):
            self.time_refreshed[k] = time_refreshed_init
            self.t_refreshed[k] = t_refreshed_init
        self.phi_values = np.zeros((self._arrays_length, 9), dtype=np.float32)


    @cy.cfunc
    def reset_size(self, n_tracked_cells: pyidx) -> cy.void:
        """
        Ensures that this can hold at least `n_tracked_cells`, resizing the internal arrays if necessary.
        Also updates `self.n_tracked_cells`.
        This can erase any data present in this object: callers must make sure this information is no longer needed.
        """
        self.n_tracked_cells = n_tracked_cells
        while self.n_tracked_cells > self._arrays_length: # dynamic resizing
            self._arrays_length *= 2
            PyMem_Free(self.ell_info)
            self.sfmin_arr = cy.cast(cy.pointer(FireBehaviorMin), PyMem_Malloc(self._arrays_length * cy.sizeof(FireBehaviorMin)))
            if not self.sfmin_arr:
                raise MemoryError()
            self.ell_info = cy.cast(cy.pointer(EllipticalInfo), PyMem_Malloc(self._arrays_length * cy.sizeof(EllipticalInfo)))
            if not self.ell_info:
                raise MemoryError()
            self.float_inputs = np.zeros((self._arrays_length, p_CellInputs), dtype=np.float32)
            PyMem_Free(self.pass1outputs)
            self.pass1outputs = cy.cast(cy.pointer(Pass1CellOutput), PyMem_Malloc(self._arrays_length * cy.sizeof(Pass1CellOutput)))
            if not self.pass1outputs:
                raise MemoryError()
            self.phi_values = np.zeros((self._arrays_length, 9), dtype=np.float32)

    def __dealloc__(self):
        PyMem_Free(self.ell_info)
        PyMem_Free(self.pass1outputs)


@cy.cfunc
@cy.exceptval(check=False)
@cy.boundscheck(False)
@cy.wraparound(False)
@cy.initializedcheck(False)
def collect_phi_values(phi: cy.float[:, :], tca: TrackedCellsArrays) -> cy.void:
    """
    Iterates over the tracked cells and stores their phi values in a cache.
    For each tracked cell, stores a row of 9 values in `tca.phi_values`,
    corresponding to the 'cross' of cells required to compute flux-limited gradients.
    Reads from `phi` and mutates `tca.phi_values`.
    """
    phi_values: cy.float[:,:] = tca.phi_values
    ell_info: cy.pointer(EllipticalInfo) = tca.ell_info
    i: pyidx
    y: pyidx
    x: pyidx
    for i in range(tca.n_tracked_cells):
        cell_index: coord_yx = ell_info[i].cell_index
        y, x = cell_index
        phi_values[i, 0] = phi[2+y, 2+x]
        phi_values[i, 1] = phi[2+y, 2+x-2]
        phi_values[i, 2] = phi[2+y, 2+x-1]
        phi_values[i, 3] = phi[2+y, 2+x+1]
        phi_values[i, 4] = phi[2+y, 2+x+2]
        phi_values[i, 5] = phi[2+y-2, 2+x]
        phi_values[i, 6] = phi[2+y-1, 2+x]
        phi_values[i, 7] = phi[2+y+1, 2+x]
        phi_values[i, 8] = phi[2+y+2, 2+x]


@cy.cfunc
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


@cy.cfunc
@cy.exceptval(check=False)
@cy.boundscheck(False)
@cy.wraparound(False)
@cy.initializedcheck(False)
def copy_tracked_cell_data(i_old: pyidx, tca_old: TrackedCellsArrays, i_new: pyidx, tca_new: TrackedCellsArrays) -> cy.void:
    # OPTIM maybe we want to use a native array directly instead of a MemoryView.
    # NOTE unrolling this loop made the code 2x faster.
    tca_new.float_inputs[i_new, 0] = tca_old.float_inputs[i_old, 0]
    tca_new.float_inputs[i_new, 1] = tca_old.float_inputs[i_old, 1]
    tca_new.float_inputs[i_new, 2] = tca_old.float_inputs[i_old, 2]
    tca_new.float_inputs[i_new, 3] = tca_old.float_inputs[i_old, 3]
    tca_new.float_inputs[i_new, 4] = tca_old.float_inputs[i_old, 4]
    tca_new.float_inputs[i_new, 5] = tca_old.float_inputs[i_old, 5]
    tca_new.float_inputs[i_new, 6] = tca_old.float_inputs[i_old, 6]
    tca_new.float_inputs[i_new, 7] = tca_old.float_inputs[i_old, 7]
    tca_new.float_inputs[i_new, 8] = tca_old.float_inputs[i_old, 8]
    tca_new.float_inputs[i_new, 9] = tca_old.float_inputs[i_old, 9]
    tca_new.float_inputs[i_new, 10] = tca_old.float_inputs[i_old, 10]
    tca_new.float_inputs[i_new, 11] = tca_old.float_inputs[i_old, 11]
    tca_new.float_inputs[i_new, 12] = tca_old.float_inputs[i_old, 12]
    tca_new.float_inputs[i_new, 13] = tca_old.float_inputs[i_old, 13]
    tca_new.float_inputs[i_new, 14] = tca_old.float_inputs[i_old, 14]
    tca_new.float_inputs[i_new, 15] = tca_old.float_inputs[i_old, 15]
    tca_new.float_inputs[i_new, 16] = tca_old.float_inputs[i_old, 16]
    tca_new.sfmin_arr[i_new] = tca_old.sfmin_arr[i_old]
    tca_new.ell_info[i_new] = tca_old.ell_info[i_old]
    # NOTE tca_old.pass1outputs does not need to be copied over given how it will get used.


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
    spot_config           : object
    refresh_frequency     : cy.float[17] # (min^-1) the frequency at which each input column needs to be refreshed.


    def __init__(self,
                 max_cells_per_timestep: float|None = 0.4,
                 buffer_width          : int|None   = 3,
                 use_wind_limit        : bool|None  = True,
                 surface_lw_ratio_model: str|None   = "behave",
                 crown_max_lw_ratio    : float|None = 1e10,
                 spot_config           : dict|None  = None,
                 inputs_refresh_freqs  : dict|None  = {}) -> cy.void:
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
            self.refresh_frequency[k] = inputs_refresh_freqs[inputs_names[k]]


@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def load_float_inputs_for_cell(
        sinputs: SpreadInputs,
        cell_index: coord_yx,
        # Where to write the data to:
        tca: TrackedCellsArrays,
        i: pyidx) -> cy.void: # NOTE maybe return the CellInputs struct instead?
    """
    Reads variables from input SpaceTimeCubes and saves them by mutating `tca.float_inputs`.
    """
    tr: pyidx[17] = tca.t_refreshed
    y, x = cell_index
    float_inputs: cy.float[:,:] = tca.float_inputs
    float_inputs[i,0] = lookup_space_time_cube_float32(sinputs.slope, (tr[0], y, x))               # rise/run
    float_inputs[i,1] = lookup_space_time_cube_float32(sinputs.aspect, (tr[1], y, x))              # degrees clockwise from North

    float_inputs[i,2] = lookup_space_time_cube_float32(sinputs.fuel_model, (tr[2], y, x))          # integer index in fm.fuel_model_table

    float_inputs[i,3] = lookup_space_time_cube_float32(sinputs.canopy_cover, (tr[3], y, x))        # 0-1
    float_inputs[i,4] = lookup_space_time_cube_float32(sinputs.canopy_height, (tr[4], y, x))       # m
    float_inputs[i,5] = lookup_space_time_cube_float32(sinputs.canopy_base_height, (tr[5], y, x))  # m
    float_inputs[i,6] = lookup_space_time_cube_float32(sinputs.canopy_bulk_density, (tr[6], y, x)) # kg/m^3

    float_inputs[i,7] = lookup_space_time_cube_float32(sinputs.wind_speed_10m, (tr[7], y, x))                # km/hr
    float_inputs[i,8] = lookup_space_time_cube_float32(sinputs.upwind_direction, (tr[8], y, x))              # degrees clockwise from North

    float_inputs[i,9] = lookup_space_time_cube_float32(sinputs.fuel_moisture_dead_1hr, (tr[9], y, x))        # kg moisture/kg ovendry weight
    float_inputs[i,10] = lookup_space_time_cube_float32(sinputs.fuel_moisture_dead_10hr, (tr[10], y, x))       # kg moisture/kg ovendry weight
    float_inputs[i,11] = lookup_space_time_cube_float32(sinputs.fuel_moisture_dead_100hr, (tr[11], y, x))      # kg moisture/kg ovendry weight
    float_inputs[i,12] = lookup_space_time_cube_float32(sinputs.fuel_moisture_live_herbaceous, (tr[12], y, x)) # kg moisture/kg ovendry weight
    float_inputs[i,13] = lookup_space_time_cube_float32(sinputs.fuel_moisture_live_woody, (tr[13], y, x))      # kg moisture/kg ovendry weight
    float_inputs[i,14] = lookup_space_time_cube_float32(sinputs.foliar_moisture, (tr[14], y, x))               # kg moisture/kg ovendry weight

    # Spread Rate Adjustments (Optional)
    float_inputs[i,15] = (lookup_space_time_cube_float32(sinputs.fuel_spread_adjustment, (tr[15], y, x))
                                 if sinputs.fuel_spread_adjustment is not None
                                 else 1.0)                                         # float >= 0.0
    float_inputs[i,16] = (lookup_space_time_cube_float32(sinputs.weather_spread_adjustment, (tr[16], y, x))
                                 if sinputs.weather_spread_adjustment is not None
                                 else 1.0)


@cy.cfunc
def list_float_inputs_cubes(sinputs: SpreadInputs) -> list[ISpaceTimeCube]:
    return [
        sinputs.slope,
        sinputs.aspect,

        sinputs.fuel_model,

        sinputs.canopy_cover,
        sinputs.canopy_height,
        sinputs.canopy_base_height,
        sinputs.canopy_bulk_density,

        sinputs.wind_speed_10m,
        sinputs.upwind_direction,

        sinputs.fuel_moisture_dead_1hr,
        sinputs.fuel_moisture_dead_10hr,
        sinputs.fuel_moisture_dead_100hr,
        sinputs.fuel_moisture_live_herbaceous,
        sinputs.fuel_moisture_live_woody,
        sinputs.foliar_moisture,

        sinputs.fuel_spread_adjustment,
        sinputs.weather_spread_adjustment,
    ]

@cy.cfunc
def default_refresh_frequency(band_duration: cy.float) -> dict:
    return {
        # Non-weather inputs default to a refresh frequency of 0.0 (never refreshed).
        'slope': 0.0,
        'aspect': 0.0,

        'fuel_model': 0.0,

        'canopy_cover': 0.0,
        'canopy_height': 0.0,
        'canopy_base_height': 0.0,
        'canopy_bulk_density': 0.0,

        # Weather inputs default to have the same refresh frequency as the base resolution of inputs.
        'wind_speed_10m': 1.0 / band_duration,
        'upwind_direction': 1.0 / band_duration,

        'fuel_moisture_dead_1hr': 1.0 / band_duration,
        'fuel_moisture_dead_10hr': 1.0 / band_duration,
        'fuel_moisture_dead_100hr': 1.0 / band_duration,
        'fuel_moisture_live_herbaceous': 1.0 / band_duration,
        'fuel_moisture_live_woody': 1.0 / band_duration,
        'foliar_moisture': 1.0 / band_duration,

        'fuel_spread_adjustment': 1.0 / band_duration,
        'weather_spread_adjustment': 1.0 / band_duration,
    }

recompute_levels_list: list = [
    100,
    100,

    100,

    100,
    100,
    100,
    100,

    10, # wind_speed_10m
    10, # upwind_direction

    100,
    100,
    100,
    100,
    100,
    100,

    100,
    100
]

@cy.cfunc
def recompute_level_for_input(input_k: pyidx) -> cy.uint:
    return recompute_levels_list[input_k]


@cy.cfunc
@cy.cdivision(True)
def refresh_inputs_if_needed(
        sinputs: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        tca: TrackedCellsArrays,
        present_time: cy.float,
        ) -> cy.uint:
    """
    Refreshes the data input columns and refresh timestamps if needed.

    Mutates `tca` and returns an integer indicating which downstream computations need to be recomputed.
    """
    recompute_level: cy.uint = 0
    stc_list: list[ISpaceTimeCube]|None = None
    k: pyidx
    for k in range(p_CellInputs):
        needs_refresh: cy.bint = (fb_opts.refresh_frequency[k] * (present_time - tca.time_refreshed[k]) > 1.0)
        if needs_refresh:
            recompute_level = max(recompute_level, recompute_level_for_input(k))
            stc_list = list_float_inputs_cubes(sinputs) if (stc_list is None) else stc_list
            space_time_cube: ISpaceTimeCube = stc_list[k]
            refresh_intvl: cy.float = 1.0 / fb_opts.refresh_frequency[k]
            time_refreshed_new: cy.float = (present_time // refresh_intvl) * refresh_intvl # NOTE REVIEW for consistency and ease of reasoning, the refresh time is always an integer multiple of the refresh interval. We might want to change this. The refresh_frequency=0 case is worth considering.
            t_refreshed_new: pyidx = int(floor(time_refreshed_new / sinputs.band_duration))
            float_inputs: cy.float[:,:] = tca.float_inputs
            y: pyidx
            x: pyidx
            i: pyidx
            for i in range(tca.n_tracked_cells):
                cell_index: coord_yx = tca.ell_info[i].cell_index
                y, x = cell_index
                float_inputs[i][k] = lookup_space_time_cube_float32(space_time_cube, (t_refreshed_new, y, x))
            tca.time_refreshed[k] = time_refreshed_new
            tca.t_refreshed[k] = t_refreshed_new
    return recompute_level



@cy.cfunc
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


@cy.cfunc
def resolve_surface_nwns_behavior(
        ci: CellInputs,
        fm_struct: FuelModel,
        ) -> FireBehaviorMin:
    """
    Computes the surface no-wind/no-slope fire behavior for a single cell.
    """
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
    # Apply fuel moisture to fuel model
    mfm: FuelModel = fm.moisturize(fm_struct, M_f)
    surface_fire_min: FireBehaviorMin = sf.calc_surface_fire_behavior_no_wind_no_slope(mfm, spread_rate_adjustment)
    return surface_fire_min



@cy.cfunc
def resolve_surface_max_behavior(
        fb_opts: FireBehaviorSettings,
        ci: CellInputs,
        fm_struct: FuelModel,
        surface_fire_min: FireBehaviorMin,
        elevation_gradient: vec_xy,
        ) -> FireBehaviorMax:
    fuel_bed_depth : cy.float = fm_struct.delta                     # ft

    # Convert from 10m wind speed to 20ft wind speed
    wind_speed_20ft: cy.float = conv.wind_speed_10m_to_wind_speed_20ft(ci.wind_speed_10m) # km/hr

    # Convert 20ft wind speed from km/hr to m/min
    wind_speed_20ft_m_min: cy.float = conv.km_hr_to_m_min(wind_speed_20ft) # m/min

    # Convert from 20ft wind speed to midflame wind speed in m/min
    midflame_wind_speed: cy.float = sf.calc_midflame_wind_speed(
                                        wind_speed_20ft_m_min,  # m/min
                                        fuel_bed_depth,         # ft
                                        conv.m_to_ft(ci.canopy_height), # ft
                                        ci.canopy_cover)           # 0-1

    # Calculate surface fire behavior in the direction of maximum spread
    surface_fire_max: FireBehaviorMax = sf.calc_surface_fire_behavior_max(
                                                surface_fire_min,
                                                midflame_wind_speed,
                                                ci.upwind_direction,
                                                ci.slope,
                                                ci.aspect,
                                                fb_opts.use_wind_limit,
                                                fb_opts.surface_lw_ratio_model)
    return surface_fire_max

@cy.cfunc
def resolve_crown_max_behavior(
        fb_opts: FireBehaviorSettings,
        ci: CellInputs,
        fm_struct: FuelModel,
        elevation_gradient: vec_xy, # NOTE not currently used but it could save some polar-to-cartesian conversion to use it.
        ) -> FireBehaviorMax:
    heat_of_combustion: cy.float = conv.Btu_lb_to_kJ_kg(fm_struct.h[0]) # kJ/kg
    estimated_fine_fuel_moisture: cy.float = ci.fuel_moisture_dead_1hr
    return cf.calc_crown_fire_behavior_max(
        ci.canopy_height,
        ci.canopy_base_height,
        ci.canopy_bulk_density,
        heat_of_combustion,
        estimated_fine_fuel_moisture,
        ci.wind_speed_10m,
        ci.upwind_direction,
        ci.slope,
        ci.aspect,
        fb_opts.crown_max_lw_ratio)


@cy.cfunc
@cy.cdivision(True)
def resolve_crowning_spread_rate_threshold(
        ci: CellInputs,
        surface_fire_max: FireBehaviorMax,
        ) -> cy.float:
    """
    Computes the surface spread rate at which crown fire occurs.
    """
    return cf.van_wagner_crowning_spread_rate_threshold(surface_fire_max, ci.canopy_base_height, ci.foliar_moisture)


@cy.cfunc
@cy.cdivision(True)
@cy.wraparound(False)
@cy.boundscheck(False)
def resolve_cell_elliptical_info(
        fb_opts: FireBehaviorSettings,
        cell_index: coord_yx,
        sinputs: SpreadInputs, # NOTE only used to call stc.get_fm_struct(fm_number). We're not really reading the rasters here.
        ci: CellInputs,
        surface_fire_min: FireBehaviorMin,
        ) -> EllipticalInfo:

    elevation_gradient: vec_xy = calc_elevation_gradient(ci.slope, ci.aspect)

    # Load the fuel model
    fm_number: pyidx = cy.cast(pyidx, ci.fuel_model_number)
    fm_struct: FuelModel = sinputs.get_fm_struct(fm_number)

    surfc_pw: PartialedEllWavelet
    crown_pw: PartialedEllWavelet
    crown_spread_rate: cy.float
    if not fm_struct.burnable:
        surfc_pw = zero_partialed_wavelet()
        crown_spread_rate = 1234.5 # arbitrary positive value - this threshold will never be reached.
        crown_pw = zero_partialed_wavelet()
    else:
        surface_fire_max: FireBehaviorMax = resolve_surface_max_behavior(fb_opts, ci, fm_struct, surface_fire_min, elevation_gradient)
        surfc_pw = pw_from_FireBehaviorMax(surface_fire_max)
        crown_spread_rate: cy.float = resolve_crowning_spread_rate_threshold(ci, surface_fire_max)
        crown_fire_max: FireBehaviorMax = resolve_crown_max_behavior(fb_opts, ci, fm_struct, elevation_gradient)
        crown_pw = pw_from_FireBehaviorMax(crown_fire_max)


    #============================================================================================
    # Build the struct
    #============================================================================================

    ell_i: EllipticalInfo = EllipticalInfo(
        cell_index = cell_index,
        slp_dz = elevation_gradient,
        surfc_wavelet = surfc_pw,
        crowning_spread_rate_threshold = crown_spread_rate,
        crown_wavelet = crown_pw
    )
    return ell_i

@cy.cfunc
@cy.initializedcheck(False)
@cy.wraparound(False)
@cy.boundscheck(False)
def refresh_caches_from_inputs_if_needed(
        sinputs: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        tca: TrackedCellsArrays,
        present_time: cy.float,
        ) -> cy.void:
    """
    If required by the refresh frequencies, refreshes inputs and recomputes the necessary downstream calcs for each tracked cell,
    such as the elliptical info and the surface no-wind/no-slope fire behavior.
    Mutates `tca`.
    """
    recompute_level: cy.uint = refresh_inputs_if_needed(sinputs, fb_opts, tca, present_time)
    i: pyidx
    float_inputs: cy.float[:,:] = tca.float_inputs
    if recompute_level >= 100:
        for i in range(tca.n_tracked_cells):
            ci: CellInputs = load_saved_CellInputs(float_inputs, i)
            fm_number: pyidx = cy.cast(pyidx, ci.fuel_model_number)
            fm_struct: FuelModel = sinputs.get_fm_struct(fm_number)
            cell_index: coord_yx = tca.ell_info[i].cell_index
            tca.sfmin_arr[i] = resolve_surface_nwns_behavior(ci, fm_struct)
    if recompute_level >= 10:
        for i in range(tca.n_tracked_cells):
            cell_index: coord_yx = tca.ell_info[i].cell_index
            ci: CellInputs = load_saved_CellInputs(float_inputs, i)
            surface_fire_min: FireBehaviorMin = tca.sfmin_arr[i]
            tca.ell_info[i] = resolve_cell_elliptical_info(fb_opts, cell_index, sinputs, ci, surface_fire_min)


@cy.cfunc
def resolve_combined_spread_behavior(
        sinputs: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        tyx: coord_tyx,
        dphi: vec_xy,
        ) -> SpreadBehavior:
    """
    Similar to resolve_cell_elliptical_info, but does a more exhaustive computation and returns the SpreadBehavior struct.
    """
    ci: CellInputs = lookup_cell_inputs(sinputs, tyx)
    elevation_gradient: vec_xy = calc_elevation_gradient(ci.slope, ci.aspect)
    dphi_st: vec_xyz = calc_phi_gradient_on_slope(dphi, elevation_gradient)
    phi_magnitude: cy.float = vu.vector_magnitude_3d(dphi_st)
    # Load the fuel model
    fm_number: pyidx = cy.cast(pyidx, ci.fuel_model_number)
    fm_struct: FuelModel = sinputs.get_fm_struct(fm_number)
    if phi_magnitude == 0.0 or not fm_struct.burnable:
        return unburned_SpreadBehavior(elevation_gradient, dphi_st)
    else:
        surface_fire_min: FireBehaviorMin = resolve_surface_nwns_behavior(ci, fm_struct)
        surface_fire_max: FireBehaviorMax = resolve_surface_max_behavior(fb_opts, ci, fm_struct, surface_fire_min, elevation_gradient)
        crown_spread_rate: cy.float = resolve_crowning_spread_rate_threshold(ci, surface_fire_max)
        sfn: SpreadBehavior = calc_fireline_normal_behavior(surface_fire_max, dphi_st)
        doesnt_crown: cy.bint = (sfn.spread_rate <= crown_spread_rate)
        if doesnt_crown:
            return sfn
        else:
            crown_fire_max: FireBehaviorMax = resolve_crown_max_behavior(fb_opts, ci, fm_struct, elevation_gradient)
            cfn: SpreadBehavior = calc_fireline_normal_behavior(crown_fire_max, dphi_st)
            combined_fire_normal: SpreadBehavior = cf.calc_combined_fire_behavior(sfn, cfn)
            return combined_fire_normal




@cy.cfunc
def load_tracked_cell_data(
        # How to resolve inputs:
        sinputs: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        cell_index: coord_yx,
        # Where to write the data to:
        tca: TrackedCellsArrays,
        i: pyidx
        ) -> cy.void:
    load_float_inputs_for_cell(sinputs, cell_index, tca, i)
    ci: CellInputs = load_saved_CellInputs(tca.float_inputs, i)
    fm_number: pyidx = cy.cast(pyidx, ci.fuel_model_number)
    fm_struct: FuelModel = sinputs.get_fm_struct(fm_number)
    surface_fire_min: FireBehaviorMin = resolve_surface_nwns_behavior(ci, fm_struct)
    tca.sfmin_arr[i] = surface_fire_min
    ell_i: EllipticalInfo = resolve_cell_elliptical_info(fb_opts, cell_index, sinputs, ci, surface_fire_min)
    tca.ell_info[i] = ell_i



@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def sync_tracked_cells_arrays(
        sinputs: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        tracked_cells: nbt.NarrowBandTracker,
        tca_old: TrackedCellsArrays,
        tca_new: TrackedCellsArrays
        ) -> cy.void:
    """
    Mutates `tca_new` so that it covers the same set of cells as `tracked_cells`,
    copying data from `tca_old` where possible,
    otherwise loading new data from `stc`.
    """
    n_tracked_cells: pyidx = tracked_cells.n_tracked_cells
    tca_new.reset_size(n_tracked_cells)
    tca_new.time_refreshed = tca_old.time_refreshed
    tca_new.t_refreshed = tca_old.t_refreshed
    cell_new: coord_yx
    i_old: pyidx = 0
    i_new: pyidx = 0
    cell_old: coord_yx = (0, 0)
    exhausted_old: cy.bint = i_old >= tca_old.n_tracked_cells
    if not(exhausted_old):
        cell_old = tca_old.ell_info[i_old].cell_index
    # This loop uses the fact that both tca_old and new_cells_itr are sorted consistently with compare_cell_indexes().
    #new_cells_itr: TrackedCellsIterator = tracked_cells_iterator(tracked_cells)
    ys_list: list = tracked_cells.ys_list
    if ys_list is not None:
        for s in ys_list:
            if s is not None:
                segm: nbt.CellsCountSegment
                for segm in s.values():
                    k: pyidx
                    y: pyidx = segm.y
                    segm_counts: cy.ushort[16] = segm.counts
                    for k in range(16):
                        if (segm_counts[k] > 0):
                            # NOTE all of the above for and if are essentially just looping over the tracked cells.
                            # This is ugly, but faster than using an Iterator pattern.
                            x: pyidx = segm.x0 + k
                            cell_new: coord_yx = (y, x)
                            while not(exhausted_old) and compare_cell_indexes(cell_old, cell_new) < 0: # cell_old is no longer tracked; just move forward.
                                i_old += 1
                                exhausted_old = i_old >= tca_old.n_tracked_cells
                                if not(exhausted_old):
                                    cell_old = tca_old.ell_info[i_old].cell_index
                            if not(exhausted_old) and (compare_cell_indexes(cell_old, cell_new) == 0): # cell_new was already tracked: copy the data.
                                copy_tracked_cell_data(i_old, tca_old, i_new, tca_new)
                            else: # cell_new was not in tca_old
                                load_tracked_cell_data(sinputs, fb_opts, cell_new, tca_new, i_new)
                            i_new += 1


@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
@cy.cdivision(True)
def runge_kutta_pass1(
        max_cells_per_timestep: cy.float,
        spatial_resolution: vec_xy,
        max_timestep: cy.float,
        tca: TrackedCellsArrays
        ) -> cy.float:
    """
    1st Runge-Kutta loop over elliptical dimensions, which:
    1. Resolves dt from the CFL condition
    2. Saves a Pass1CellOutput struct for each cell.

    Returns the resolved `dt` and mutates `tca.pass1outputs`.
    Reads only `tca.ell_info`.
    """
    ell_info: cy.pointer[EllipticalInfo] = tca.ell_info
    pass1outputs: cy.pointer[Pass1CellOutput] = tca.pass1outputs
    (dx, dy) = spatial_resolution
    # The following will be useful to compute dt based on the CFL constraint.
    # It is more convenient to first compute dt_inv, the reciprocal of dt;
    # dt_inv = 0 represents an infinite dt. We will later enforce that dt <= max_timestep.
    dt_inv: cy.float = 0.0
    C_dx: cy.float = max_cells_per_timestep * dx
    C_dy: cy.float = max_cells_per_timestep * dy
    # Now looping over tracked cells:
    phi_values: cy.float[:,:] = tca.phi_values
    dx_inv: cy.float = 1.0 / dx
    dy_inv: cy.float = 1.0 / dy
    for i in range(tca.n_tracked_cells):
        ell_i: EllipticalInfo = ell_info[i]
        cell_index: coord_yx = ell_i.cell_index
        # y: pyidx = cell_index[0]
        # x: pyidx = cell_index[1]
        # dphi: vec_xy = calc_phi_gradient_approx(phi, dx, dy, x, y)
        dphi_dx: cy.float = (phi_values[i, 3] - phi_values[i, 2]) * dx_inv / 2.0
        dphi_dy: cy.float = (phi_values[i, 7] - phi_values[i, 6]) * dy_inv / 2.0
        dphi: vec_xy = (dphi_dx, dphi_dy)
        dphi_norm2: cy.float = (dphi_dx * dphi_dx) + (dphi_dy * dphi_dy)
        dphi_dt_flim: cy.float
        if dphi_norm2 > 0: # Most common case.
            dphi_dx_flim: cy.float = calc_dphi_flim_x(phi_values[i, 0], phi_values[i, 1], phi_values[i, 2], phi_values[i, 3], phi_values[i, 4]) * dx_inv
            dphi_dy_flim: cy.float = calc_dphi_flim_y(phi_values[i, 0], phi_values[i, 5], phi_values[i, 6], phi_values[i, 7], phi_values[i, 8]) * dy_inv
            dphi_dt: cy.float = dphi_dt_from_elliptical(ell_i, dphi)
            dphi_dt_correction: cy.float = (dphi_dx * dphi_dx_flim + dphi_dy * dphi_dy_flim) / dphi_norm2
            dphi_dt_flim = (dphi_dt * dphi_dt_correction)
            # Checking the CFL condition and updating dt_inv if needed (which will be rare).
            # The code is written in this way to be fast, but it's not trivial that it's correct; proof below.
            # The CFL constraint is defined as the requirement that |Ux*dt| <= C*dx and |Uy*dt| <= C*dy,
            # in which U := (Ux, Uy) is the front-normal spread rate vector in the horizontal plane,
            # and C := max_cells_per_timestep.
            # Recall that we could express U as follows: U: vec_xy = scale_2d(-dphi_dt/dphi_norm2, dphi),
            # which follows from the facts that dphi_dt = - dot2d(U, dphi) and that U is by definition positively proportional to dphi.
            # In particular, Ux = - dphi_dx * dphi_dt / dphi_norm2.
            # Our contraint (from Ux) thus becomes:
            # |dt* dphi_dx * dphi_dt / dphi_norm2| <= C * dx
            # Recalling that dt_inv := 1/dt and rearranging yields:
            # dt_inv * (C * dx) * dphi_norm2 >= |dphi_dx * dphi_dt|
            if (dt_inv * (dphi_norm2 * C_dx) < abs(dphi_dt * dphi_dx)): # dt is too large given Ux.
                dt_inv = abs(dphi_dt * dphi_dx) / (dphi_norm2 * C_dx)
            # And similarly for Uy:
            if (dt_inv * (dphi_norm2 * C_dy) < abs(dphi_dt * dphi_dy)): # dt is too large given Uy.
                dt_inv = abs(dphi_dt * dphi_dy) / (dphi_norm2 * C_dy)
        else:
            dphi_dt_flim = 0.0
        pass1outputs[i] = Pass1CellOutput( # INTRO
            cell_index = cell_index,
            dphi = dphi,
            dphi_dt_flim = dphi_dt_flim,
            phi_old = phi_values[i, 0]
            )
    dt_inv = max(dt_inv, 1.0/max_timestep) # (dt <= max_timestep) iff (dt_inv >= 1/max_timestep).
    dt: cy.float = 1.0/dt_inv
    return dt


@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def update_phi_star(
        tca: TrackedCellsArrays,
        dt: cy.float,
        phs: cy.float[:,:]
        ) -> cy.void:
    """
    Mutates the phs ('phi star') matrix, by using the dt and dphi/dt computed in the 1st Runge-Kutta pass.
    To be called between Runge-Kutta passes.
    """
    pass1outputs: cy.pointer[Pass1CellOutput] = tca.pass1outputs
    i: pyidx
    for i in range(tca.n_tracked_cells):
        pass1out: Pass1CellOutput = pass1outputs[i]
        cell_index: coord_yx = pass1out.cell_index
        y: pyidx = cell_index[0]
        x: pyidx = cell_index[1]
        dphi_dt_0i: cy.float = pass1out.dphi_dt_flim
        phs[2+y, 2+x] = pass1out.phi_old + (dt * dphi_dt_0i)


@cy.cclass
class BurnedCellInfo: # Using an Extension Type instead of a struct because it's convenient to store in Python data structures like lists and dicts.
    """
    This data structure simply records information about a burned cell.
    """
    cell_index: coord_yx
    toa: cy.float # Time Of Arrival
    dphi: vec_xy # Gradient of Phi field, indicative of front direction.
    from_spotting: cy.bint # Whether spotting is what caused the cell to ignite.


@cy.cfunc
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


@cy.cfunc
@cy.cdivision(True)
@cy.wraparound(False)
@cy.boundscheck(False)
def runge_kutta_pass2(
        fb_opts: FireBehaviorSettings,
        spatial_resolution: vec_xy,
        start_time: cy.float,
        dt: cy.float,
        tca: TrackedCellsArrays,
        phi: cy.float[:, :]
        ) -> list[BurnedCellInfo]:
    """
    2nd Runge-Kutta loop, which:
    1. Updates the phi matrix
    2. Identifies cells that have just burned and returns them in a list.
    Reads from `tca` and `phs`, and mutates `phi`.
    """
    (dx, dy) = spatial_resolution
    n_tracked_cells: pyidx = tca.n_tracked_cells
    ell: cy.pointer[EllipticalInfo] = tca.ell_info
    pass1outputs: cy.pointer[Pass1CellOutput] = tca.pass1outputs
    spread_burned_cells: list[BurnedCellInfo] = []
    i: pyidx
    phs_values: cy.float[:,:] = tca.phi_values
    dx_inv: cy.float = 1.0 / dx
    dy_inv: cy.float = 1.0 / dy
    for i in range(n_tracked_cells):
        ell_i: EllipticalInfo = ell[i]
        cell_index: coord_yx = ell_i.cell_index
        y: pyidx = cell_index[0]
        x: pyidx = cell_index[1]
        #dphi: vec_xy = calc_phi_gradient_approx(phs, dx, dy, x, y)
        dphi_dx: cy.float = (phs_values[i, 3] - phs_values[i, 2]) * dx_inv / 2
        dphi_dy: cy.float = (phs_values[i, 7] - phs_values[i, 6]) * dy_inv / 2
        dphi: vec_xy = (dphi_dx, dphi_dy)
        dphi_norm2: cy.float = (dphi_dx * dphi_dx) + (dphi_dy * dphi_dy)
        dphi_dt: cy.float
        dphi_dt_1i: cy.float
        if dphi_norm2 > 0: # Most common case.
            dphi_dx_flim: cy.float = calc_dphi_flim_x(phs_values[i, 0], phs_values[i, 1], phs_values[i, 2], phs_values[i, 3], phs_values[i, 4]) * dx_inv
            dphi_dy_flim: cy.float = calc_dphi_flim_y(phs_values[i, 0], phs_values[i, 5], phs_values[i, 6], phs_values[i, 7], phs_values[i, 8]) * dy_inv
            dphi_dt_correction: cy.float = (dphi_dx * dphi_dx_flim + dphi_dy * dphi_dy_flim) / dphi_norm2
            dphi_dt = dphi_dt_from_elliptical(ell_i, dphi)
            dphi_dt_1i = (dphi_dt * dphi_dt_correction)
        else:
            dphi_dt_1i = 0.0
        dphi_dt_0i: cy.float = pass1outputs[i].dphi_dt_flim
        phi_old: cy.float = pass1outputs[i].phi_old
        phi_new: cy.float = phi_old + 0.5 * dt * (dphi_dt_0i + dphi_dt_1i)
        phi[2+y, 2+x] = phi_new
        i_just_burned: cy.bint = (phi_old * phi_new) < 0 # Phi can only ever decrease, therefore if these are of opposite signs, the cell has just burned.
        if i_just_burned:
            spread_burned_cells.append(new_BurnedCellInfo(
                cell_index,
                toa = start_time + dt * phi_old / (phi_old - phi_new),
                dphi = pass1outputs[i].dphi, # Using the 1st-pass phi gradient. FIXME to be consistent with the toa, we might want to average this with the dphi from the 2nd pass.
                from_spotting=False))
    return spread_burned_cells


@cy.cfunc
@cy.cdivision(True)
@cy.wraparound(False)
@cy.boundscheck(False)
def process_spread_burned_cells(
        sinputs: SpreadInputs,
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
        t: pyidx = int(toa // sinputs.band_duration)
        cell_index: coord_yx = burned_cell.cell_index
        y, x = cell_index
        tyx: coord_tyx = (t, y, x)
        # Re-compute the spread behavior. It's OK to re-compute it because a cell burning is a relatively rare event.
        fb: SpreadBehavior = resolve_combined_spread_behavior(sinputs, fb_opts, tyx, burned_cell.dphi)
        # Write to outputs
        fire_type_matrix[y,x]          = fb.fire_type
        spread_rate_matrix[y,x]        = fb.spread_rate
        spread_direction_matrix[y,x]   = vu.spread_direction_vector_to_angle(fb.spread_direction)
        fireline_intensity_matrix[y,x] = fb.fireline_intensity
        flame_length_matrix[y,x]       = fb.flame_length
        time_of_arrival_matrix[y,x]    = toa

        # Cast firebrands and update spot_ignitions
        if fb_opts.spot_config:
            spot_from_burned_cell(
                fb_opts.spot_config,
                sinputs,
                fire_type_matrix,
                random_generator,
                y,
                x,
                fb,
                toa,
                spot_ignitions)


@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def reset_phi_star(
        tca: TrackedCellsArrays,
        spot_ignited: list[BurnedCellInfo],
        phi_star_matrix: cy.float[:,:],
        phi_matrix: cy.float[:,:]
    ) -> cy.void:
    """
    Efficiently updates `phi_star_matrix` to match `phi_matrix`,
    by copying only the values of cells where phi have changed.
    Mutates `phi_star_matrix`, reading from `tca.pass1outputs`, `spot_ignited` and `phi_matrix`.
    """
    y: pyidx
    x: pyidx
    i: pyidx
    for i in range(tca.n_tracked_cells):
        y, x = tca.pass1outputs[i].cell_index
        phi_star_matrix[2+y, 2+x] = phi_matrix[2+y, 2+x]
    tcb: BurnedCellInfo
    for tcb in spot_ignited:
        y, x = tcb.cell_index
        phi_star_matrix[2+y, 2+x] = phi_matrix[2+y, 2+x]


@cy.cfunc
@cy.wraparound(False)
@cy.boundscheck(False)
def ignite_from_spotting(spot_ignitions: SortedDict, output_matrices, stop_time: cy.float) -> list[BurnedCellInfo]:
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
                if phi_matrix[2+y, 2+x] > 0.0: # Not burned yet
                    phi_matrix[2+y, 2+x]             = -1.0
                    time_of_arrival_matrix[y,x] = ignition_time # FIXME: REVIEW Should I use stop_time instead?
                    ignited.append(new_BurnedCellInfo(
                        cell_index,
                        toa = ignition_time,
                        dphi = (0.0, 0.0),
                        from_spotting = True
                    ))
                    # FIXME: I need to calculate and store the fire_behavior values for these cells
    return ignited

@cy.cfunc
@cy.exceptval(check=False)
@cy.inline
def encode_cell_index(y: pyidx, x: pyidx) -> object:
    """
    Encodes a (y, x) tuple into a single Python integer object.
    This enables a more efficient memory layout than a tuple of Python integers.
    """
    return ((cy.cast(cy.ulong, y) << 32) + cy.cast(cy.ulong, x))

@cy.cfunc
@cy.exceptval(check=False)
@cy.inline
def decode_cell_index(cell_index: object) -> coord_yx:
    ci: cy.ulong = cell_index
    y: pyidx = ci >> 32
    x: pyidx = cy.cast(cy.uint, ci) # NOTE I would prefer writing this as (ci & 0xFFFFFFFF), but Cython compilation makes it slower.
    return (y, x)

@cy.cfunc
def route_cell_to_diff(
        frontier_cells_old: set,
        phi: cy.float[:, :],
        frontier_additions: set,
        frontier_removals: set,
        y: pyidx,
        x: pyidx,
        ) -> cy.void:
    """
    Determines whether the given cell was just added or removed from the frontier cells,
    mutating the sets `frontier_additions` and `frontier_removals` accordingly.
    Idempotent.
    """
    set_elem: object = encode_cell_index(y, x)
    if is_frontier_cell(phi, y, x):
        if not (set_elem in frontier_cells_old):
            frontier_additions.add(set_elem)
    else:
        if (set_elem in frontier_cells_old):
            frontier_removals.add(set_elem)


@cy.cfunc
def diff_frontier_cells(
        frontier_cells_old: set,
        phi: cy.float[:, :],
        spread_ignited: list[BurnedCellInfo],
        spot_ignited: list[BurnedCellInfo],
        rows: pyidx,
        cols: pyidx
    ) -> object:
    """
    Computes the diff between the old frontier cells and the new frontier cells, based on newly burned cells.
    Returns a (cells_added, cells_dropped) tuple of sets containing cell indices encoded by `encode_cell_index`.
    """
    frontier_additions: set = set()
    frontier_removals: set = set()
    l: list[BurnedCellInfo]
    for l in [spread_ignited, spot_ignited]: # NOTE we accept two lists instead of one to avoid paying the cost of concatenating them.
        ci: BurnedCellInfo
        y: pyidx
        x: pyidx
        for ci in l:
            y, x = ci.cell_index
            # NOTE only in the neighborhood of a burned cell can there be changes to frontier cells membership.
            route_cell_to_diff(frontier_cells_old, phi, frontier_additions, frontier_removals, y, x)
            route_cell_to_diff(frontier_cells_old, phi, frontier_additions, frontier_removals, y-1, x)
            route_cell_to_diff(frontier_cells_old, phi, frontier_additions, frontier_removals, y+1, x)
            route_cell_to_diff(frontier_cells_old, phi, frontier_additions, frontier_removals, y, x-1)
            route_cell_to_diff(frontier_cells_old, phi, frontier_additions, frontier_removals, y, x+1)
    return (frontier_additions, frontier_removals)


@cy.cfunc
def apply_frontier_diff(frontier_cells_old: set, frontier_additions: set, frontier_removals: set) -> set:
    frontier_cells_new: set = frontier_cells_old.copy()
    for set_elem in frontier_additions:
        frontier_cells_new.add(set_elem)
    for set_elem in frontier_removals:
        frontier_cells_new.discard(set_elem)
    return frontier_cells_new


@cy.cfunc
def update_tracked_cells_with_frontier_diff(
        tracked_cells: object,
        frontier_cells_added: set,
        frontier_cells_dropped: set,
        buffer_width: pyidx
    ) -> object:
    """
    TODO: Add docstring
    """
    # Determine which frontier cells have been added or dropped
    # OPTIM these set differences are 75% of the runtime of this function.
    # We might be able to resolve additions and drops in a faster way by iterating over burned cells.
    cell                  : object
    # Increment reference counters for all cells within buffer_width of the added frontier cells
    y: pyidx
    x: pyidx
    for cell in frontier_cells_added:
        y, x = decode_cell_index(cell)
        nbt.inc_square_around(tracked_cells, y, x, buffer_width)
    # Decrement reference counters for all cells within buffer_width of the dropped frontier cells
    for cell in frontier_cells_dropped:
        y, x = decode_cell_index(cell)
        nbt.dec_square_around(tracked_cells, y, x, buffer_width)
    # Return updated tracked cells
    return tracked_cells


@cy.cfunc
@cy.cdivision(True)
def spread_one_timestep(
        sim_state: dict,
        sinputs: SpreadInputs,
        fb_opts: FireBehaviorSettings,
        max_timestep: cy.float,
    ) -> dict:
    """
    Spreads the fire for one iteration using the level-set method, returning an updated `sim_state`.

    """
    output_matrices: dict = sim_state["output_matrices"]
    spot_ignitions: object = sim_state["spot_ignitions"]
    random_generator: BufferedRandGen = sim_state["random_generator"]
    phi: cy.float[:,:] = output_matrices["phi"]
    phs: cy.float[:,:] = output_matrices["phi_star"]

    # Insert missing tracked cells.
    tracked_cells: nbt.NarrowBandTracker = sim_state["tracked_cells"]
    tca_old: TrackedCellsArrays = sim_state["_tracked_cells_arrays_old"]
    tca: TrackedCellsArrays = sim_state["_tracked_cells_arrays"]
    start_time: cy.float = sim_state["simulation_time"]
    sync_tracked_cells_arrays(sinputs, fb_opts, tracked_cells, tca_old, tca)

    refresh_caches_from_inputs_if_needed(sinputs, fb_opts, tca, start_time)

    collect_phi_values(phi, tca)
    dt = runge_kutta_pass1(fb_opts.max_cells_per_timestep, sinputs.spatial_resolution, max_timestep, tca)

    # Now that dt is known, update phi_star_matrix.
    update_phi_star(tca, dt, phs)
    stop_time: cy.float = start_time + dt

    collect_phi_values(phs, tca)
    spread_burned_cells: list[BurnedCellInfo] = runge_kutta_pass2(fb_opts, sinputs.spatial_resolution, start_time, dt, tca, phi)

    # Side-effects of the burned cells (outputs etc.).
    process_spread_burned_cells(sinputs, fb_opts, output_matrices, spot_ignitions, random_generator, spread_burned_cells)
    # TODO REVIEW It is a questionable choice to call this function AFTER process_spread_burned_cells;
    # it may be more sensible to ignite the spotting cells first and then to process them all.
    spot_ignited: list[BurnedCellInfo] = ignite_from_spotting(spot_ignitions, output_matrices, stop_time)

    # Save the new phi_matrix values in phi_star_matrix
    reset_phi_star(tca, spot_ignited, phs, phi)

    # Update the sets of frontier cells and tracked cells based on the updated phi matrix
    frontier_cells_old: set  = sim_state["frontier_cells"] # OPTIM maybe it would be more efficient to use a binary array here.
    frontier_additions, frontier_removals = diff_frontier_cells(
        frontier_cells_old,
        phi,
        spread_burned_cells,
        spot_ignited,
        sinputs.rows,
        sinputs.cols
    )
    frontier_cells_new: set = apply_frontier_diff(frontier_cells_old, frontier_additions, frontier_removals)
    tracked_cells_new : object = update_tracked_cells_with_frontier_diff(tracked_cells, frontier_additions, frontier_removals, fb_opts.buffer_width)

    # Return the updated world state
    sim_state_new = {
        "simulation_time" : stop_time,
        "output_matrices" : output_matrices,
        "frontier_cells"  : frontier_cells_new,
        "tracked_cells"   : tracked_cells_new,
        # Purposefully swapping the tracked_cells_arrays
        "_tracked_cells_arrays": tca_old,
        "_tracked_cells_arrays_old": tca,

        "spot_ignitions"  : spot_ignitions,
        "random_generator": random_generator,
    }
    return sim_state_new



def spread_fire_with_phi_field(space_time_cubes: dict, output_matrices: dict, cube_resolution: tuple,
                               start_time: float, max_duration: float|None = None,
                               max_cells_per_timestep: float|None = 0.4, buffer_width: int|None = 3,
                               use_wind_limit: bool|None = True, surface_lw_ratio_model: str|None = "behave",
                               crown_max_lw_ratio: float|None = 1e10, spot_ignitions: dict|None = {},
                               spot_config: dict|None = None, inputs_refresh_freqs: dict|None = {}):
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
    - inputs_refresh_freqs         :: dictionary from input name to refresh frequency in 1/min (Optional).
                                      0 means never refresh. Weather inputs default to 1/band_duration,
                                      whereas non-weather inputs default to 0.

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
    - spot_ignitions       :: dictionary of (ignition_time -> ignited_cells) (only included when spotting is used)
    - random_generator     :: numpy.random.Generator object (only included when spotting is used)
    """
    # Define the provided, required, and optional keys for space_time_cubes
    provided_cubes = set(space_time_cubes.keys())
    required_cubes = {
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
    optional_cubes = {
        "fuel_spread_adjustment",
        "weather_spread_adjustment",
    }

    # Ensure that all required_cubes are present in provided_cubes
    if not provided_cubes.issuperset(required_cubes):
        raise ValueError("The space_time_cubes dictionary is missing these required keys: "
                         + str(required_cubes.difference(provided_cubes)))

    # Ensure that only required_cubes and optional_cubes are present in provided_cubes
    if not (required_cubes | optional_cubes).issuperset(provided_cubes):
        raise ValueError("The space_time_cubes dictionary contains these unused keys: "
                         + str(provided_cubes.difference((required_cubes | optional_cubes))))

    # Ensure that all space_time_cubes values are ISpaceTimeCube objects
    for cube in space_time_cubes.values():
        if not(isinstance(cube, ISpaceTimeCube)):
            raise ValueError("All values in the space_time_cubes dictionary must be SpaceTimeCube or "
                             + "LazySpaceTimeCube objects. See pyretechnics.space_time_cube for more information.")

    # Define the provided, required, and optional keys for output_matrices
    provided_matrices = set(output_matrices.keys())
    required_matrices = {
        "phi",
        "fire_type",
        "spread_rate",
        "spread_direction",
        "fireline_intensity",
        "flame_length",
        "time_of_arrival",
    }
    optional_matrices = set() # NOTE: Leaving this here in case we add some later

    # Ensure that all required_matrices are present in output_matrices
    if not provided_matrices.issuperset(required_matrices):
        raise ValueError("The output_matrices dictionary is missing these required keys: "
                         + str(required_matrices.difference(provided_matrices)))

    # Ensure that only required_matrices and optional_matrices are present in provided_matrices
    if not (required_matrices | optional_matrices).issuperset(provided_matrices):
        raise ValueError("The output_matrices dictionary contains these unused keys: "
                         + str(provided_matrices.difference((required_matrices | optional_matrices))))

    # Ensure that all output_matrices values are 2D Numpy arrays
    for matrix in output_matrices.values():
        if not(isinstance(matrix, np.ndarray) and np.ndim(matrix) == 2):
            raise ValueError("All values in the output_matrices dictionary must be 2D Numpy arrays.")

    # Extract simulation dimensions
    (bands, rows, cols) = space_time_cubes["slope"].shape

    # Extract simulation resolution
    (band_duration, cell_height, cell_width) = cube_resolution

    # Ensure that all space_time_cubes have the same spatial resolution
    for cube in space_time_cubes.values():
        if cube.shape != (bands, rows, cols):
            raise ValueError("The space_time_cubes must all share the same spatial resolution.")

    sinputs: SpreadInputs = make_SpreadInputs(cube_resolution, space_time_cubes)
    fb_opts: FireBehaviorSettings = FireBehaviorSettings(
        max_cells_per_timestep = max_cells_per_timestep,
        buffer_width           = buffer_width,
        use_wind_limit         = use_wind_limit,
        surface_lw_ratio_model = surface_lw_ratio_model,
        crown_max_lw_ratio     = crown_max_lw_ratio,
        spot_config            = spot_config,
        inputs_refresh_freqs   = {**default_refresh_frequency(band_duration), **inputs_refresh_freqs},
    )

    # Ensure that space_time_cubes and output_matrices have the same spatial resolution
    for (label, matrix) in output_matrices.items():
        if label == 'phi':
            if matrix.shape != (rows + 4, cols + 4):
                raise ValueError("The phi matrix must be padded by 2 rows and 2 columns compared to the simulation dimensions.")
        else:
            if matrix.shape != (rows, cols):
                raise ValueError("The space_time_cubes and output_matrices must share the same spatial resolution.")

    # Ensure that all cube resolution values are positive
    if band_duration <= 0.0 or cell_height <= 0.0 or cell_width <= 0.0:
        raise ValueError("The cube_resolution tuple may only contain positive values.")

    # Calculate the cube duration
    cube_duration = bands * band_duration

    # Ensure that start_time exists within the temporal bounds of the space_time_cubes
    if not(0.0 <= start_time < cube_duration):
        raise ValueError("The start_time falls outside of the temporal bounds of the space_time_cubes.")

    # Calculate the max stop time
    max_stop_time = start_time + max_duration if max_duration else cube_duration

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
    spot_igns = SortedDict(spot_ignitions)

    start_t: pyidx = int(start_time // band_duration)
    sim_state = {
        "simulation_time" : start_time,
        "output_matrices" : output_matrices,
        "frontier_cells"  : frontier_cells,
        "tracked_cells"   : tracked_cells,
        # Purposefully swapping the tracked_cells_arrays
        "_tracked_cells_arrays": TrackedCellsArrays(start_time, start_t), # It's OK not to be in sync - spread_one_timestep will solve this.
        "_tracked_cells_arrays_old": TrackedCellsArrays(start_time, start_t),
        "spot_ignitions"  : spot_igns,
        "random_generator": random_generator,
    }

    # Spread the fire until an exit condition is reached
    # FIXME: I don't think the "no burnable cells" condition can ever be met currently.
    simulation_time: cy.float = sim_state["simulation_time"]

    while(simulation_time < max_stop_time and (nbt.nonempty_tracked_cells(tracked_cells) or len(spot_ignitions) > 0)):
        # Compute max_timestep based on the remaining time in the temporal band and simulation
        remaining_time_in_band       = band_duration - simulation_time % band_duration
        remaining_time_in_simulation = max_stop_time - simulation_time
        max_timestep                 = min(remaining_time_in_band, remaining_time_in_simulation)
        # Spread fire one timestep
        sim_state = spread_one_timestep(sim_state, sinputs, fb_opts, max_timestep)

        # Reset spread inputs
        simulation_time  = sim_state["simulation_time"]

        tracked_cells    = sim_state["tracked_cells"]
        spot_ignitions   = sim_state["spot_ignitions"]

    output_matrices  = sim_state["output_matrices"]
    # Remove the temporary copy of the phi matrix from output_matrices # NOTE REVIEW I think phi_star would be better off in sim_state.
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
