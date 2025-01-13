# [[file:../../org/pyretechnics.org::van-wagner-critical-fireline-intensity][van-wagner-critical-fireline-intensity]]
# cython: profile=True
import cython
if cython.compiled:
    from cython.cimports.pyretechnics.math import sqrt
else:
    from math import sqrt
import cython as cy


@cy.ccall
@cy.profile(False)
def van_wagner_critical_fireline_intensity(canopy_base_height: cy.float, foliar_moisture: cy.float) -> cy.float:
    """
    Returns the critical fireline intensity (kW/m) given:
    - canopy_base_height :: m
    - foliar_moisture    :: kg moisture/kg ovendry weight

    Constants used:
    460.0 = heat-of-ignition :: kJ/kg
    0.01 = empirical estimate for C in Van Wagner 1977 (eq. 4)
    """
    H = 460.0 + 2600.0 * foliar_moisture
    v: cy.float = 0.01 * canopy_base_height * H
    return v * sqrt(v) # NOTE the original implementation did pow(v, 1.5). Refactored to use sqrt because it's faster.
# van-wagner-critical-fireline-intensity ends here
# [[file:../../org/pyretechnics.org::van-wagner-crown-fire-initiation][van-wagner-crown-fire-initiation]]
import cython
if cython.compiled:
    import cython.cimports.pyretechnics.cy_types as py_types
else:
    import pyretechnics.py_types as py_types
import cython as cy

@cy.ccall
@cy.cdivision(True)
@cy.exceptval(check=False)
def van_wagner_crowning_spread_rate_threshold(
        surface_fire_max: py_types.FireBehaviorMax,
        canopy_base_height: cy.float,
        foliar_moisture: cy.float,
        ) -> cy.float:
    """
    Computes the surface spread rate above which crown fire occurs.
    """
    surfc_max_fli: cy.float = surface_fire_max.max_fireline_intensity
    crown_fli: cy.float = van_wagner_critical_fireline_intensity(canopy_base_height, foliar_moisture)
    crowning_spread_rate_threshold: cy.float = (surface_fire_max.max_spread_rate * crown_fli / surfc_max_fli # NOTE this uses the fact that fireline intensity is proportional to spread rate.
                                                if surfc_max_fli > 0
                                                else 0)
    return crowning_spread_rate_threshold

@cy.ccall
def van_wagner_crown_fire_initiation(
        surface_fireline_intensity: cy.float, 
        canopy_cover: cy.float, 
        canopy_base_height: cy.float, 
        foliar_moisture: cy.float
        ) -> cy.bint:
    """
    Returns True if the surface fire transitions to a crown fire or False otherwise given:
    - surface_fireline_intensity :: kW/m
    - canopy_cover               :: 0-1
    - canopy_base_height         :: m
    - foliar_moisture            :: kg moisture/kg ovendry weight
    """
    return (
        surface_fireline_intensity > 0.0
        and
        canopy_cover > 0.4
        and
        surface_fireline_intensity >= van_wagner_critical_fireline_intensity(canopy_base_height, foliar_moisture)
    )
# van-wagner-crown-fire-initiation ends here
# [[file:../../org/pyretechnics.org::cruz-active-crown-fire-spread-rate][cruz-active-crown-fire-spread-rate]]
if cython.compiled:
    from cython.cimports.pyretechnics.math import exp, pow
else:
    from math import exp, pow


@cy.ccall
@cy.exceptval(check=False)
def cruz_active_crown_fire_spread_rate(wind_speed_10m: cy.float, canopy_bulk_density: cy.float, estimated_fine_fuel_moisture: cy.float) -> cy.float:
    """
    Returns the active crown fire spread rate (m/min) given:
    - wind_speed_10m                                   :: km/hr
    - canopy_bulk_density                              :: kg/m^3
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") :: kg moisture/kg ovendry weight
    """
    return (11.02
            * pow(wind_speed_10m, 0.90)
            * pow(canopy_bulk_density, 0.19)
            * exp(-17.0 * estimated_fine_fuel_moisture))
# cruz-active-crown-fire-spread-rate ends here
# [[file:../../org/pyretechnics.org::van-wagner-critical-spread-rate][van-wagner-critical-spread-rate]]
@cy.ccall
@cy.cdivision(True)
@cy.inline
def van_wagner_critical_spread_rate(canopy_bulk_density: cy.float) -> cy.float:
    """
    Returns the critical spread rate (m/min) given:
    - canopy_bulk_density :: kg/m^3
    """
    return 3.0 / canopy_bulk_density
# van-wagner-critical-spread-rate ends here
# [[file:../../org/pyretechnics.org::cruz-passive-crown-fire-spread-rate][cruz-passive-crown-fire-spread-rate]]
if cython.compiled:
    from cython.cimports.pyretechnics.math import exp
else:
    from math import exp

@cy.ccall
@cy.cdivision(True)
@cy.exceptval(-1)
def cruz_passive_crown_fire_spread_rate(active_spread_rate: cy.float, critical_spread_rate: cy.float) -> cy.float:
    """
    Returns the passive crown fire spread rate (m/min) given:
    - active_spread_rate   :: m/min
    - critical_spread_rate :: m/min
    """
    return active_spread_rate * exp(-active_spread_rate / critical_spread_rate)
# cruz-passive-crown-fire-spread-rate ends here
# [[file:../../org/pyretechnics.org::cruz-crown-fire-spread-info][cruz-crown-fire-spread-info]]
CrownSpreadInfo = cy.struct(
    fire_type = cy.int,
    spread_rate = cy.float,
    critical_spread_rate = cy.float # INTRO the spread rate separating passive and active crowning. 
)

@cy.ccall
def cruz_crown_fire_spread_info(
        wind_speed_10m: cy.float,
        canopy_bulk_density: cy.float,
        estimated_fine_fuel_moisture: cy.float
        ) -> CrownSpreadInfo:
    """
    Given these inputs:
    - wind_speed_10m                                   :: km/hr
    - canopy_bulk_density                              :: kg/m^3
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") :: kg moisture/kg ovendry weight

    return a dictionary containing these keys:
    - fire_type            :: "passive_crown" or "active_crown"
    - spread_rate          :: m/min
    - critical_spread_rate :: m/min
    """
    active_spread_rate   = cruz_active_crown_fire_spread_rate(wind_speed_10m,
                                                              canopy_bulk_density,
                                                              estimated_fine_fuel_moisture) # m/min
    critical_spread_rate = van_wagner_critical_spread_rate(canopy_bulk_density) # m/min
    if (active_spread_rate > critical_spread_rate):
        return CrownSpreadInfo(
            fire_type = 3,  # FIXME NAMED CONSTANT
            spread_rate = active_spread_rate,
            critical_spread_rate = critical_spread_rate
        )
    else:
        return CrownSpreadInfo(
            fire_type = 2,  # FIXME NAMED CONSTANT
            spread_rate = cruz_passive_crown_fire_spread_rate(active_spread_rate, critical_spread_rate),
            critical_spread_rate = critical_spread_rate
        )
# cruz-crown-fire-spread-info ends here
# [[file:../../org/pyretechnics.org::crown-fireline-intensity][crown-fireline-intensity]]
# NOTE: heat_of_combustion is h from the fuel models (generally 8000 Btu/lb)
# NOTE: ELMFIRE hard-codes heat_of_combustion to 18000 kJ/kg = 7738.6 Btu/lb
@cy.ccall
@cy.exceptval(check=False)
@cy.inline
def calc_crown_fireline_intensity(
        crown_spread_rate: cy.float, 
        canopy_bulk_density: cy.float,
        canopy_height: cy.float,
        canopy_base_height: cy.float, 
        heat_of_combustion: cy.float
        ) -> cy.float:
    """
    Returns the crown fireline intensity (Btu/ft/s OR kW/m) given:
    - crown_spread_rate                                             :: ft/min  OR m/min
    - canopy_bulk_density                                           :: lb/ft^3 OR kg/m^3
    - canopy_height                                                 :: ft      OR m
    - canopy_base_height                                            :: ft      OR m
    - heat_of_combustion                                            :: Btu/lb  OR kJ/kg

    (ft/min * lb/ft^3 * ft * Btu/lb)/60 = (Btu/ft/min)/60 = Btu/ft/s
    OR
    (m/min * kg/m^3 * m * kJ/kg)/60 = (kJ/m*min)/60 = kJ/m*s = kW/m
    """
    canopy_height_difference: cy.float = canopy_height - canopy_base_height
    return (crown_spread_rate * canopy_bulk_density * canopy_height_difference * heat_of_combustion) / 60.0
# crown-fireline-intensity ends here
# [[file:../../org/pyretechnics.org::crown-fire-eccentricity][crown-fire-eccentricity]]
if cython.compiled:
    from cython.cimports.pyretechnics.math import sqrt
else:
    from math import sqrt
import pyretechnics.conversion as conv

# Parameters for the linear model that computes LoW from wind speed.
LoW_intercept = cy.declare(cy.float, 1.0) 
LoW_slope_per_km_hr = cy.declare(cy.float, 0.07767140120267868) 
assert LoW_slope_per_km_hr == conv.km_hr_to_mph(0.125) # The original formula used wind speeds in mph 


@cy.ccall
@cy.exceptval(check=False)
@cy.inline
def crown_length_to_width_ratio(wind_speed_10m: cy.float, max_length_to_width_ratio: cy.float=1e10) -> cy.float:
    """
    Calculate the length_to_width_ratio of the crown fire front using eq. 9 from
    Rothermel 1991 given:
    - wind_speed_10m            :: km/hr (aligned with the slope-tangential plane)
    - max_length_to_width_ratio :: float > 0.0 (Optional)
    """
    # FIXME why is wind speed not in SI units???
    length_to_width_ratio: cy.float = LoW_intercept + LoW_slope_per_km_hr * wind_speed_10m
    return min(length_to_width_ratio, max_length_to_width_ratio)


@cy.ccall
@cy.cdivision(True)
def crown_fire_eccentricity(length_to_width_ratio: cy.float) -> cy.float:
    """
    Calculate the eccentricity (E) of the crown fire front using eq. 8 from
    Albini and Chase 1980 given:
    - L/W :: (1: circular spread, > 1: elliptical spread)
    """
    LoW2: cy.float = length_to_width_ratio * length_to_width_ratio
    return sqrt(LoW2 - 1.0) / length_to_width_ratio
# crown-fire-eccentricity ends here
# [[file:../../org/pyretechnics.org::crown-fire-behavior-max][crown-fire-behavior-max]]
import cython
if cython.compiled:
    import cython.cimports.pyretechnics.conversion as conv
    import cython.cimports.pyretechnics.cy_types as py_types
    import cython.cimports.pyretechnics.surface_fire1 as sf
    import cython.cimports.pyretechnics.vector_utils as vu
else:
    import pyretechnics.conversion as conv
    import pyretechnics.py_types as py_types
    import pyretechnics.surface_fire1 as sf
    import pyretechnics.vector_utils as vu
import numpy as np
import cython as cy

@cy.profile(False)
@cy.ccall
def calc_crown_fire_behavior_max(
        canopy_height: cy.float, 
        canopy_base_height: cy.float, 
        canopy_bulk_density: cy.float, 
        heat_of_combustion: cy.float,
        estimated_fine_fuel_moisture: cy.float, 
        wind_speed_10m: cy.float, 
        upwind_direction: cy.float,
        slope: cy.float, 
        aspect: cy.float, 
        crown_max_lw_ratio: cy.float=1e10,
        ) -> py_types.FireBehaviorMax:
    """
    Given these inputs:
    - canopy_height                                    :: m
    - canopy_base_height                               :: m
    - canopy_bulk_density                              :: kg/m^3
    - heat_of_combustion                               :: kJ/kg
    - estimated_fine_fuel_moisture (M_f[0] "dead-1hr") :: kg moisture/kg ovendry weight
    - wind_speed_10m                                   :: km/hr
    - upwind_direction                                 :: degrees clockwise from North
    - slope                                            :: rise/run
    - aspect                                           :: degrees clockwise from North
    - crown_max_lw_ratio                               :: float > 0.0 (Optional)

    return a dictionary containing these keys:
    - max_fire_type          :: "passive_crown" or "active_crown"
    - max_spread_rate        :: m/min
    - max_spread_direction   :: (x, y, z) unit vector
    - max_fireline_intensity :: kW/m
    - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
    - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
    - critical_spread_rate   :: m/min
    """
    # Reverse the provided wind and slope directions
    downwind_direction = conv.opposite_direction(upwind_direction)
    upslope_direction  = conv.opposite_direction(aspect)
    # Project wind and slope vectors onto the slope-tangential plane
    # FIXME let's just have these vectors as arguments to the function instead of re-computing them.
    vectors: sf.ProjectedVectors      = sf.project_wind_and_slope_vectors_3d(wind_speed_10m, downwind_direction, slope, upslope_direction)
    wind_vector_3d: py_types.vec_xyz  = vectors.wind_vector_3d  # km/hr
    slope_vector_3d: py_types.vec_xyz = vectors.slope_vector_3d # rise/run
    # Determine the max spread direction
    wind_speed_10m_3d: cy.float = vu.vector_magnitude_3d(wind_vector_3d)      # km/hr
    max_spread_direction: py_types.vec_xyz
    if wind_speed_10m_3d > 0.0:
        max_spread_direction = vu.as_unit_vector_3d(wind_vector_3d)       # unit vector in the 3D downwind direction
    elif slope > 0.0:
        max_spread_direction = vu.as_unit_vector_3d(slope_vector_3d) # unit vector in the 3D upslope direction
    else:
        max_spread_direction = (0.0,1.0,0.0) # default: North
    # Calculate the crown fire behavior in the max spread direction
    spread_info: CrownSpreadInfo    = cruz_crown_fire_spread_info(wind_speed_10m_3d, canopy_bulk_density, estimated_fine_fuel_moisture)
    spread_rate: cy.float           = spread_info.spread_rate # m/min
    fireline_intensity: cy.float    = calc_crown_fireline_intensity(spread_rate, canopy_bulk_density, canopy_height,
                                                          canopy_base_height, heat_of_combustion) # kW/m
    length_to_width_ratio: cy.float = crown_length_to_width_ratio(wind_speed_10m_3d, crown_max_lw_ratio) # unitless
    eccentricity: cy.float          = crown_fire_eccentricity(length_to_width_ratio) # unitless
    return py_types.FireBehaviorMax(
        spread_info.fire_type,
        spread_rate,
        max_spread_direction, # unit vector
        fireline_intensity,
        -1.0, # NOTE max_flame_length is not provided, as in the original unoptimized code.
        length_to_width_ratio,
        eccentricity,
        spread_info.critical_spread_rate,
    )
# crown-fire-behavior-max ends here
# [[file:../../org/pyretechnics.org::crown-fire-behavior-in-direction][crown-fire-behavior-in-direction]]
import numpy as np


def calc_crown_fire_behavior_in_direction(crown_fire_max, spread_direction):
    """
    Given these inputs:
    - crown_fire_max     :: dictionary of max crown fire behavior values
      - max_fire_type          :: "passive_crown" or "active_crown"
      - max_spread_rate        :: m/min
      - max_spread_direction   :: (x, y, z) unit vector
      - max_fireline_intensity :: kW/m
      - length_to_width_ratio  :: unitless (1: circular spread, > 1: elliptical spread)
      - eccentricity           :: unitless (0: circular spread, > 0: elliptical spread)
      - critical_spread_rate   :: m/min
    - spread_direction   :: 3D unit vector on the slope-tangential plane

    return a dictionary containing these keys:
    - fire_type          :: "passive_crown" or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    """
    # Unpack max crown fire behavior values
    max_fire_type          = crown_fire_max["max_fire_type"]
    max_spread_rate        = crown_fire_max["max_spread_rate"]
    max_spread_direction   = crown_fire_max["max_spread_direction"]
    max_fireline_intensity = crown_fire_max["max_fireline_intensity"]
    eccentricity           = crown_fire_max["eccentricity"]
    critical_spread_rate   = crown_fire_max["critical_spread_rate"]
    # Calculate cos(w), where w is the offset angle between these unit vectors on the slope-tangential plane
    cos_w = np.dot(max_spread_direction, spread_direction)
    # Calculate adjustment due to the offset angle from the max spread direction
    adjustment = (1.0 - eccentricity) / (1.0 - eccentricity * cos_w)
    # Adjust the spread rate (possibly switching from an active to passive crown fire)
    spread_rate = max_spread_rate * adjustment
    if spread_rate > critical_spread_rate:
        # Max spread rate was active and directional spread rate remains active
        return {
            "fire_type"         : "active_crown",
            "spread_rate"       : spread_rate,
            "spread_direction"  : np.asarray(spread_direction),
            "fireline_intensity": max_fireline_intensity * adjustment,
        }
    elif max_fire_type == "passive_crown":
        # Max spread rate was passive and directional spread rate remains passive
        return {
            "fire_type"         : "passive_crown",
            "spread_rate"       : spread_rate,
            "spread_direction"  : np.asarray(spread_direction),
            "fireline_intensity": max_fireline_intensity * adjustment,
        }
    else:
        # Max spread rate was active and directional spread rate has become passive
        return {
            "fire_type"         : "passive_crown",
            "spread_rate"       : cruz_passive_crown_fire_spread_rate(spread_rate, critical_spread_rate),
            "spread_direction"  : np.asarray(spread_direction),
            "fireline_intensity": max_fireline_intensity * adjustment,
        }
# crown-fire-behavior-in-direction ends here
# [[file:../../org/pyretechnics.org::combined-fire-behavior][combined-fire-behavior]]


@cy.ccall
def calc_combined_fire_behavior(surface_fire_behavior: py_types.SpreadBehavior, crown_fire_behavior: py_types.SpreadBehavior) -> py_types.SpreadBehavior:
    """
    Given these inputs:
    - surface_fire_behavior :: dictionary of surface fire behavior values
      - fire_type              :: "surface"
      - spread_rate            :: m/min
      - spread_direction       :: (x, y, z) unit vector
      - fireline_intensity     :: kW/m
      - flame_length           :: m
    - crown_fire_behavior   :: dictionary of crown fire behavior values
      - fire_type              :: "passive_crown" or "active_crown"
      - spread_rate            :: m/min
      - spread_direction       :: (x, y, z) unit vector
      - fireline_intensity     :: kW/m

    return a dictionary containing these keys:
    - fire_type          :: "surface", "passive_crown", or "active_crown"
    - spread_rate        :: m/min
    - spread_direction   :: (x, y, z) unit vector
    - fireline_intensity :: kW/m
    - flame_length       :: m
    """
    # Unpack the surface fire behavior values
    surface_spread_rate        = surface_fire_behavior.spread_rate        # m/min
    surface_spread_direction: py_types.vec_xyz = surface_fire_behavior.spread_direction   # (x, y, z) unit vector
    surface_fireline_intensity = surface_fire_behavior.fireline_intensity # kW/m
    # Unpack the crown fire behavior values
    crown_fire_type          = crown_fire_behavior.fire_type          # "passive_crown" or "active_crown"
    crown_spread_rate        = crown_fire_behavior.spread_rate        # m/min
    crown_spread_direction: py_types.vec_xyz = crown_fire_behavior.spread_direction   # (x, y, z) unit vector
    crown_fireline_intensity = crown_fire_behavior.fireline_intensity # kW/m
    
    surfc_dphi_dt = surface_fire_behavior.dphi_dt
    crown_dphi_dt = crown_fire_behavior.dphi_dt
    dphi_dt: cy.float
    if surfc_dphi_dt < crown_dphi_dt:
        dphi_dt = surfc_dphi_dt
    else:
        dphi_dt = crown_dphi_dt

    # Determine whether the surface or crown fire has the fastest spread rate
    if surface_spread_rate == 0.0:
        # Independent crown fire (This is probably user error.)
        return crown_fire_behavior
    elif crown_spread_rate == 0.0:
        # No crown fire
        return surface_fire_behavior
    elif surface_spread_rate > crown_spread_rate:
        # Surface fire spreads faster
        combined_fireline_intensity = (surface_fireline_intensity
                                       + crown_fireline_intensity * surface_spread_rate / crown_spread_rate)
        return py_types.SpreadBehavior(
            dphi_dt,
            crown_fire_type,
            surface_spread_rate,
            surface_spread_direction,
            combined_fireline_intensity,
            sf.calc_flame_length(combined_fireline_intensity),
        )
    else:
        # Crown fire spreads faster
        combined_fireline_intensity = (surface_fireline_intensity * crown_spread_rate / surface_spread_rate
                                       + crown_fireline_intensity)
        return py_types.SpreadBehavior(
            dphi_dt,
            crown_fire_type,
            crown_spread_rate,
            crown_spread_direction,
            combined_fireline_intensity,
            sf.calc_flame_length(combined_fireline_intensity),
        )
# combined-fire-behavior ends here
