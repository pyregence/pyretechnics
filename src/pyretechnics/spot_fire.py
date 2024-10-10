# [[file:../../org/pyretechnics.org::expected-ember-production][expected-ember-production]]
from math import sqrt


def expected_ember_production(fire_behavior, elevation_gradient, cube_resolution, ember_production_rate=0.03):
    """
    Return the expected number of embers produced by an entire cell when it burns given:
    - fire_behavior         :: dictionary of surface or crown fire behavior values
      - fire_type                :: "unburned", "surface", "passive_crown", or "active_crown"
      - spread_rate              :: m/min
      - spread_direction         :: (x, y, z) unit vector on the slope-tangential plane
      - fireline_intensity       :: kW/m
      - flame_length             :: m
    - elevation_gradient    :: tuple with these fields:
      - dz_dx                    :: rise/run
      - dz_dy                    :: rise/run
    - cube_resolution       :: tuple with these fields:
      - band_duration            :: minutes
      - cell_height              :: meters
      - cell_width               :: meters
    - ember_production_rate :: embers/kJ
    """
    if fire_behavior["spread_rate"] == 0.0:
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

        (dz_dx, dz_dy) = elevation_gradient                      # (rise/run, rise/run)
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
