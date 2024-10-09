# [[file:../../org/pyretechnics.org::expected-ember-production][expected-ember-production]]
from math import sqrt
from pyretechnics.eulerian_level_set import calc_elevation_gradient


def expected_ember_production(space_time_cubes, space_time_coordinate, cube_resolution, fire_behavior,
                              ember_production_rate=0.03):
    """
    Return the expected number of embers produced at the specified space_time_coordinate given:
    - space_time_cubes          :: dictionary of (Lazy)SpaceTimeCube objects with these cell types
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
    - space_time_coordinate     :: (t,y,x)
    - cube_resolution           :: tuple with these fields:
      - band_duration                 :: minutes
      - cell_height                   :: meters
      - cell_width                    :: meters
    - fire_behavior             :: dictionary of surface or crown fire behavior values
      - fire_type                     :: "unburned", "surface", "passive_crown", or "active_crown"
      - spread_rate                   :: m/min
      - spread_direction              :: (x, y, z) unit vector on the slope-tangential plane
      - fireline_intensity            :: kW/m
      - flame_length                  :: m
    - ember_production_rate     :: embers/kJ
    """
    if fire_behavior["fire_type"] == "unburned":
        return 0.0
    else:
        #================================================================================================
        # Calculate the heat output per unit area
        #================================================================================================

        spread_rate          = fire_behavior["spread_rate"]            # m/min
        fireline_intensity   = fire_behavior["fireline_intensity"]     # kW/m
        heat_output_per_area = 60.0 * fireline_intensity / spread_rate # kJ/m^2 (reaction_intensity * residence_time)

        #================================================================================================
        # Calculate the slope-adjusted cell area
        #================================================================================================

        (t, y, x)      = space_time_coordinate
        slope          = space_time_cubes["slope"].get(t, y, x)  # rise/run
        aspect         = space_time_cubes["aspect"].get(t, y, x) # degrees clockwise from North
        (dz_dx, dz_dy) = calc_elevation_gradient(slope, aspect)  # (rise/run, rise/run)
        slope_factor   = sqrt(1.0 + dz_dx ** 2.0 + dz_dy ** 2.0) # unitless
        cell_height    = cube_resolution[1]                      # meters
        cell_width     = cube_resolution[2]                      # meters
        cell_area      = cell_height * cell_width * slope_factor # m^2

        #================================================================================================
        # Calculate the expected number of embers produced in this cell
        #================================================================================================

        cell_heat_output = heat_output_per_area * cell_area         # kJ
        ember_count      = cell_heat_output * ember_production_rate # number of embers
        return ember_count
# expected-ember-production ends here
