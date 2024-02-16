#==============================================================
# Fire Spread Functions
#==============================================================

# FIXME: stub
def burn_cells(layer_dict, ignited_cells_set, outputs_list):
    """
    1. Load data for the ignited area without a buffer
    2. Burn all the ignited cells
    3. Return a dictionary of the output strings to their corresponding 2D matrices
    """
    return None

# FIXME: stub
def spread(layer_dict, ignited_cells_set, outputs_list, stop_condition):
    """
    1. Load data for the ignited area plus a buffer size (user-specified?)
    2. Perform constant spread out over the landscape in all directions
       - Run surface, crown, and spot equations per ignited cell
    3. Record the time-of-arrival (ignition-time) in each cell as it spreads
    4. Load more data whenever the buffer extent is exceeded or stop spreading if no more data is available
    5. Continue until a stop condition is met
    6. Return a dictionary of the output strings to their corresponding 2D matrices (start with time-of-arrival)
    """
    return None
