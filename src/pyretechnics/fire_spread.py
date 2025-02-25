# [[file:../../org/pyretechnics.org::fire-spread-functions][fire-spread-functions]]
# TODO: Compare this with numpy.random (is there a generator object that I can use here?)
import random
import numpy as np


def supply_layer(layer_dict, layer_name, layer_shape, layer_type, initial_value, overwrite_layers):
    layer = layer_dict.get(layer_name)
    if layer:
        if overwrite_layers:
            return layer
        else:
            return layer.copy()
    else:
        return np.full(layer_shape, initial_value, dtype=layer_type)


# FIXME: stub
def burn_cells(world_state):
    """
    1. Validate inputs
       - Make sure that all required parameters are present in config_dict
       - Make sure that all required entries are present in layer_dict
       - Check for entries in output_layer_dict; if present, copy/use as new output layers
         - Use time of arrival raster to identify the burn scar(s); if missing, we are simulating point ignitions
       - Make sure that ignited_cells_set is a set of integer 3-tuples
    2. Iterate through all cells in (rows, cols), compute the maximum spread rate and severity values, and store them in output arrays
       - Run surface and crown equations per cell
    3. Return a dictionary of the outputs
    """
    simulation_clock  = world_state["simulation_clock"]
    config_dict       = world_state["config_dict"]
    input_layer_dict  = world_state["input_layer_dict"]
    output_layer_dict = world_state.get("output_layer_dict", {})
    ignited_cells_set = world_state["ignited_cells_set"]

    # The random seed affects input layer perturbations.
    if config_dict.get("random_seed"):
        random.seed(config_dict.get("random_seed"))

    return None


# FIXME: stub
def grow_fire(world_state):
    """
    1. Validate inputs
       - Make sure that all required parameters are present in config_dict
       - Make sure that all required entries are present in layer_dict
       - Check for entries in output_layer_dict; if present, copy/use as new output layers
         - Use time of arrival raster to identify the burn scar(s); if missing, we are simulating point ignitions
       - Make sure that ignited_cells_set is a set of integer 3-tuples
    2. Iterate through all cells in ignited_cells_set, compute the perimeter-oriented spread rate and severity values, and store them in output arrays
    3. Perform constant spread out over the landscape in all directions
       - Run surface, crown, and spot equations per ignited cell
    4. Record the time_of_arrival in each cell as it spreads
    5. Continue until a stop condition is met (e.g., max_burn_duration or max_burned_area)
    6. Return a dictionary of the outputs
    """
    simulation_clock  = world_state["simulation_clock"]
    config_dict       = world_state["config_dict"]
    input_layer_dict  = world_state["input_layer_dict"]
    output_layer_dict = world_state.get("output_layer_dict", {})
    ignited_cells_set = world_state["ignited_cells_set"]

    # The random seed affects input layer perturbations and spotting calculations.
    if config_dict.get("random_seed"):
        random.seed(config_dict.get("random_seed"))

    # GridFire initializes the following 2D arrays for its spread algorithm:
    #
    #   :burn-time-matrix                [float32] time_of_arrival (with -1 for values in the interior of the burn scar)
    #   :eccentricity-matrix             [float32] -1 in burn scar, 0 otherwise
    #   :fireline-intensity-matrix       [float32] -1 in burn scar, 0 otherwise
    #   :fire-spread-matrix              [float32] +1 in burn scar, 0 otherwise
    #   :fire-type-matrix                [float32] -1 in burn scar, 0 otherwise
    #   :firebrand-count-matrix          [ int32 ]  0 everywhere (when spotting params are passed)
    #   :flame-length-matrix             [float32] -1 in burn scar, 0 otherwise
    #   :directional-flame-length-matrix [float32] -1 in burn scar, 0 otherwise (when compute-directional-values? = true)
    #   :max-spread-direction-matrix     [float32] -1 in burn scar, 0 otherwise
    #   :max-spread-rate-matrix          [float32] -1 in burn scar, 0 otherwise
    #   :modified-time-matrix            [ int32 ]  0 everywhere
    #   :residence-time-matrix           [float32] -1 in burn scar, 0 otherwise (when compute-directional-values? = true)
    #   :reaction-intensity-matrix       [float32] -1 in burn scar, 0 otherwise (when compute-directional-values? = true)
    #   :spot-matrix                     [float32]  0 everywhere (when spotting params are passed)
    #   :spread-rate-matrix              [float32] -1 in burn scar, 0 otherwise
    #   :spread-rate-sum-matrix          [float32]  0 everywhere (when compute-directional-values? = true)
    #   :travel-lines-matrix             [ int16 ]  0 everywhere
    #   :x-magnitude-sum-matrix          [float32]  0 everywhere (when compute-directional-values? = true)
    #   :y-magnitude-sum-matrix          [float32]  0 everywhere (when compute-directional-values? = true)
    #
    # Equivalent Pyretechnics 2D arrays in output_layer_dict:
    #
    #   ========================= Output 2D Arrays =========================
    #   time_of_arrival                  :burn-time-matrix
    #   max_surface_spread_direction     :max-spread-direction-matrix
    #   max_crown_spread_direction       :max-spread-direction-matrix
    #   perimeter_spread_direction       N/A
    #   max_surface_spread_rate          :max-spread-rate-matrix
    #   max_crown_spread_rate            :max-spread-rate-matrix
    #   perimeter_spread_rate            :spread-rate-matrix
    #   max_surface_fireline_intensity   :fireline-intensity-matrix
    #   max_crown_fireline_intensity     :fireline-intensity-matrix
    #   perimeter_fireline_intensity     :fireline-intensity-matrix
    #   max_surface_flame_length         :flame-length-matrix
    #   max_crown_flame_length           :flame-length-matrix
    #   perimeter_flame_length           :directional-flame-length-matrix
    #   fire_type                        :fire-type-matrix
    #   firebrand_ignition               :spot-matrix
    #
    #   ======================== Internal 2D Arrays ========================
    #   surface_eccentricity             :eccentricity-matrix
    #   crown_eccentricity               :eccentricity-matrix
    #
    # GridFire 2D arrays that are not needed by Pyretechnics:
    #
    #   ======================== Replaced 2D Arrays ========================
    #   :fire-spread-matrix              time_of_arrival
    #   :residence-time-matrix           max_surface_fireline_intensity, max_crown_fireline_intensity, surface_eccentricity, crown_eccentricity
    #   :reaction-intensity-matrix       max_surface_fireline_intensity, max_crown_fireline_intensity, surface_eccentricity, crown_eccentricity
    #
    #   ======================== Internal 2D Arrays ========================
    #   :modified-time-matrix            \
    #   :travel-lines-matrix              |
    #   :spread-rate-sum-matrix           |-- for its 2D spread algorithm
    #   :x-magnitude-sum-matrix           |
    #   :y-magnitude-sum-matrix          /

    (num_timesteps, num_rows, num_cols) = config_dict["simulation_shape"]
    layer_shape = (num_rows, num_cols)
    overwrite_outputs = config_dict["overwrite_outputs"]

    output_layer_dict = {
        "eulerian_level_set_phi_field"  : supply_layer(output_layer_dict, "eulerian_level_set_phi_field"  , layer_shape, "float16", np.nan, overwrite_outputs),
        "time_of_arrival"               : supply_layer(output_layer_dict, "time_of_arrival"               , layer_shape, "float32", np.nan, overwrite_outputs),
        "max_surface_spread_direction"  : supply_layer(output_layer_dict, "max_surface_spread_direction"  , layer_shape, "float16", np.nan, overwrite_outputs),
        "max_crown_spread_direction"    : supply_layer(output_layer_dict, "max_crown_spread_direction"    , layer_shape, "float16", np.nan, overwrite_outputs),
        "perimeter_spread_direction"    : supply_layer(output_layer_dict, "perimeter_spread_direction"    , layer_shape, "float16", np.nan, overwrite_outputs),
        "max_surface_spread_rate"       : supply_layer(output_layer_dict, "max_surface_spread_rate"       , layer_shape, "float16", np.nan, overwrite_outputs),
        "max_crown_spread_rate"         : supply_layer(output_layer_dict, "max_crown_spread_rate"         , layer_shape, "float16", np.nan, overwrite_outputs),
        "perimeter_spread_rate"         : supply_layer(output_layer_dict, "perimeter_spread_rate"         , layer_shape, "float16", np.nan, overwrite_outputs),
        "surface_eccentricity"          : supply_layer(output_layer_dict, "surface_eccentricity"          , layer_shape, "float16", np.nan, overwrite_outputs),
        "crown_eccentricity"            : supply_layer(output_layer_dict, "crown_eccentricity"            , layer_shape, "float16", np.nan, overwrite_outputs),
        "max_surface_fireline_intensity": supply_layer(output_layer_dict, "max_surface_fireline_intensity", layer_shape, "float32", np.nan, overwrite_outputs),
        "max_crown_fireline_intensity"  : supply_layer(output_layer_dict, "max_crown_fireline_intensity"  , layer_shape, "float32", np.nan, overwrite_outputs),
        "perimeter_fireline_intensity"  : supply_layer(output_layer_dict, "perimeter_fireline_intensity"  , layer_shape, "float32", np.nan, overwrite_outputs),
        "max_surface_flame_length"      : supply_layer(output_layer_dict, "max_surface_flame_length"      , layer_shape, "float16", np.nan, overwrite_outputs),
        "max_crown_flame_length"        : supply_layer(output_layer_dict, "max_crown_flame_length"        , layer_shape, "float16", np.nan, overwrite_outputs),
        "perimeter_flame_length"        : supply_layer(output_layer_dict, "perimeter_flame_length"        , layer_shape, "float16", np.nan, overwrite_outputs),
        "fire_type"                     : supply_layer(output_layer_dict, "fire_type"                     , layer_shape, "uint8"  ,      0, overwrite_outputs),
        "firebrand_ignition"            : supply_layer(output_layer_dict, "firebrand_ignition"                , layer_shape, "bool8"  ,  False, overwrite_outputs),
    }

    # RESUME at [[file:~/code/sig-gis/gridfire/src/gridfire/fire_spread.clj::(defn- run-loop]]
    # TODO: Investigate ELMFIRE's inputs to determine if we are missing anything needed by its API.

    return None
# fire-spread-functions ends here
