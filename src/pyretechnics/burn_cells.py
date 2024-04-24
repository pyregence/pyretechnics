# [[file:../../org/Pyretechnics.org::burn-cells][burn-cells]]
from math import sqrt, atan2, degrees
from fuel_models import fuel_models_precomputed, moisturize
import surface_fire as sf
from conversion import wind_speed_10m_to_wind_speed_20ft, m_to_ft, mps_to_fpm

# NOTE: Name change (no more !)
# NOTE: Make inputs["fuel_spread_adjustment"] and inputs["weather_spread_adjustment"] default to a function that always returns 1.0
# NOTE: The wind layers are now using cartesian rather than polar coordinates
def compute_max_in_situ_values(inputs, outputs, t, i, j):
    """
    Computes and saves fire-spread behavior values for the space-time cell at location (t, i, j).
    """
    # 2D Inputs
    slope                         = inputs["slope"](i, j)
    aspect                        = inputs["aspect"](i, j)
    fuel_model_number             = inputs["fuel_model"](i, j)
    canopy_cover                  = inputs["canopy_cover"](i, j)
    canopy_height                 = inputs["canopy_height"](i, j)
    canopy_base_height            = inputs["canopy_base_height"](i, j)
    canopy_bulk_density           = inputs["canopy_bulk_density"](i, j)
    fuel_spread_adjustment        = inputs["fuel_spread_adjustment"](i, j)
    # 3D Inputs
    wind_speed_10m_x              = inputs["wind_speed_10m_x"](t, i, j)
    wind_speed_10m_y              = inputs["wind_speed_10m_y"](t, i, j)
    fuel_moisture_dead_1hr        = inputs["fuel_moisture_dead_1hr"](t, i, j)
    fuel_moisture_dead_10hr       = inputs["fuel_moisture_dead_10hr"](t, i, j)
    fuel_moisture_dead_100hr      = inputs["fuel_moisture_dead_100hr"](t, i, j)
    fuel_moisture_live_herbaceous = inputs["fuel_moisture_live_herbaceous"](t, i, j)
    fuel_moisture_live_woody      = inputs["fuel_moisture_live_woody"](t, i, j)
    foliar_moisture               = inputs["foliar_moisture"](t, i, j)
    weather_spread_adjustment     = inputs["weather_spread_adjustment"](t, i, j)
    # 2D Output Arrays (FIXME: match these to output_layer_dict)
    max_spread_rate_matrix        = outputs["max_spread_rate"]
    max_spread_direction_matrix   = outputs["max_spread_direction"]
    spread_rate_matrix            = outputs["spread_rate"]
    flame_length_matrix           = outputs["flame_length"]
    fire_line_intensity_matrix    = outputs["fire_line_intensity"]
    fire_type_matrix              = outputs["fire_type"]
    modified_time_matrix          = outputs["modified_time"]
    eccentricity_matrix           = outputs["eccentricity"]
    residence_time_matrix         = outputs["residence_time"]
    reaction_intensity_matrix     = outputs["reaction_intensity"]
    # Surface Fire Min (TODO: memoize rothermel_surface_fire_spread_no_wind_no_slope)
    fuel_model                    = fuel_models_precomputed[fuel_model_number]
    fuel_moisture                 = [fuel_moisture_dead_1hr,
                                     fuel_moisture_dead_10hr,
                                     fuel_moisture_dead_100hr,
                                     0.0, # fuel_moisture_dead_herbaceous
                                     fuel_moisture_live_herbaceous,
                                     fuel_moisture_live_woody]
    moisturized_fuel_model        = moisturize(fuel_model, fuel_moisture)
    surface_fire_min              = sf.rothermel_surface_fire_spread_no_wind_no_slope(moisturized_fuel_model)
    # Midflame Wind Speed
    wind_speed_10m                = sqrt(wind_speed_10m_x ** 2.0 + wind_speed_10m_y ** 2.0)
    wind_speed_20ft               = wind_speed_10m_to_wind_speed_20ft(wind_speed_10m)
    wind_from_direction           = (90.0 - degrees(atan2(wind_speed_10m_y, wind_speed_10m_x)) % 360.0) % 360.0
    wind_adj_factor               = sf.wind_adjustment_factor(fuel_model["delta"], m_to_ft(canopy_height), canopy_cover)
    midflame_wind_speed           = mps_to_fpm(wind_speed_20ft * wind_adj_factor)
    # Surface Fire Max
    spread_rate_adjustment        = fuel_spread_adjustment * weather_spread_adjustment
    surface_fire_max              = sf.rothermel_surface_fire_spread_max(surface_fire_min,
                                                                         midflame_wind_speed,
                                                                         wind_from_direction,
                                                                         slope,
                                                                         aspect,
                                                                         spread_rate_adjustment)
    # Max Surface Intensity
    max_spread_rate               = surface_fire_max["max_spread_rate"]
    residence_time                = surface_fire_min["residence_time"]
    reaction_intensity            = surface_fire_min["reaction_intensity"]
    # max_spread_direction          = surface_fire_max["max_spread_direction"]
    # eccentricity                  = surface_fire_max["eccentricity"]
    #=======================================================================================
    # NOTE: The calculations to determine the fireline_normal_spread_rate have been elided.
    #       Consider ending this function here and making another function to compute the
    #       surface/crown values based on the provided perimeter spread direction.
    #=======================================================================================
    max_flame_depth               = anderson_flame_depth(max_spread_rate, residence_time)
    max_surface_intensity         = byram_fire_line_intensity(reaction_intensity, max_flame_depth)
    # ?
    if van_wagner_crown_fire_initiation(canopy_cover,
                                        canopy_base_height,
                                        foliar_moisture,
                                        max_surface_intensity):
        crown_spread_max        = cruz_crown_fire_spread(wind_speed_20ft,
                                                         canopy_bulk_density,
                                                         fuel_moisture_dead_1hr)
        crown_type              = 2.0 if (crown_spread_max < 0.0) else 3.0 # 2=passive, 3=active
        crown_spread_max        = abs(crown_spread_max)
        crown_eccentricity      = eccentricity if (max_spread_rate > crown_spread_max) else crown_fire_eccentricity(wind_speed_20ft)
        #=======================================================================================
        # NOTE: The calculations to determine the fireline_normal_spread_rate have been elided.
        #=======================================================================================
        max_crown_intensity     = crown_fire_line_intensity(crown_spread_max,
                                                            canopy_bulk_density,
                                                            (canopy_height - canopy_base_height),
                                                            fuel_model["h"][0])
        tot_fire_line_intensity = surface_intensity + crown_intensity
        max_spread_rate         = max(max_spread_rate, crown_spread_max)

        # LET BODY: RESUME HERE
        (t/mset! max_spread_rate_matrix i j max_spread_rate)
        (t/mset! max_spread_direction_matrix i j max_spread_direction)
        (t/mset! eccentricity_matrix i j crown_eccentricity)
        (t/mset! modified_time_matrix i j (inc t))
        # compute_directional_values?                      (:compute_directional_values? inputs)
        (when compute_directional_values?
          (t/mset! residence_time_matrix i j residence_time)
          (t/mset! reaction_intensity_matrix i j reaction_intensity))
        (store_if_max! spread_rate_matrix i j crown_fln_spread_rate)
        (store_if_max! flame_length_matrix i j (byram_flame_length tot_fire_line_intensity))
        (store_if_max! fire_line_intensity_matrix i j tot_fire_line_intensity)
        (store_if_max! fire_type_matrix i j crown_type))
      (do
        (t/mset! max_spread_rate_matrix i j max_spread_rate)
        (t/mset! max_spread_direction_matrix i j max_spread_direction)
        (t/mset! eccentricity_matrix i j eccentricity)
        (t/mset! modified_time_matrix i j (inc t))
        (when compute_directional_values?
          (t/mset! residence_time_matrix i j residence_time)
          (t/mset! reaction_intensity_matrix i j reaction_intensity))
        (store_if_max! spread_rate_matrix i j fireline_normal_spread_rate)
        (store_if_max! flame_length_matrix i j (byram_flame_length surface_intensity))
        (store_if_max! fire_line_intensity_matrix i j surface_intensity)
        (store_if_max! fire_type_matrix i j 1.0)))))
# burn-cells ends here
