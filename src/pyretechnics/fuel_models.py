# [[file:../../org/pyretechnics.org::fuel-model-compact-table][fuel-model-compact-table]]
import cython as cy

# Lookup table including entries for each of the Anderson 13 and Scott & Burgan 40 fuel models.
#
# The fields have the following meanings:
#   {fuel-model-number : [name, delta, M_x-dead, h, w_o, sigma]}
#
# where:
#   w_o   = [  w_o-dead-1hr,   w_o-dead-10hr,   w_o-dead-100hr,   w_o-live-herbaceous,   w_o-live-woody]
#   sigma = [sigma-dead-1hr, sigma-dead-10hr, sigma-dead-100hr, sigma-live-herbaceous, sigma-live-woody]
fuel_model_compact_table = {
    # Anderson 13:
    # Grass and Grass-dominated (short-grass,timber-grass-and-understory,tall-grass)
    1   : ["R01", 1.0, 12, 8, [0.0340, 0.0000, 0.0000, 0.0000, 0.0000], [3500.0,   0.0,  0.0,    0.0,    0.0]],
    2   : ["R02", 1.0, 15, 8, [0.0920, 0.0460, 0.0230, 0.0230, 0.0000], [3000.0, 109.0, 30.0, 1500.0,    0.0]],
    3   : ["R03", 2.5, 25, 8, [0.1380, 0.0000, 0.0000, 0.0000, 0.0000], [1500.0,   0.0,  0.0,    0.0,    0.0]],
    # Chaparral and Shrubfields (chaparral,brush,dormant-brush-hardwood-slash,southern-rough)
    4   : ["R04", 6.0, 20, 8, [0.2300, 0.1840, 0.0920, 0.2300, 0.0000], [2000.0, 109.0, 30.0, 1500.0,    0.0]],
    5   : ["R05", 2.0, 20, 8, [0.0460, 0.0230, 0.0000, 0.0920, 0.0000], [2000.0, 109.0,  0.0, 1500.0,    0.0]],
    6   : ["R06", 2.5, 25, 8, [0.0690, 0.1150, 0.0920, 0.0000, 0.0000], [1750.0, 109.0, 30.0,    0.0,    0.0]],
    7   : ["R07", 2.5, 40, 8, [0.0520, 0.0860, 0.0690, 0.0170, 0.0000], [1750.0, 109.0, 30.0, 1550.0,    0.0]],
    # Timber Litter (closed-timber-litter,hardwood-litter,timber-litter-and-understory)
    8   : ["R08", 0.2, 30, 8, [0.0690, 0.0460, 0.1150, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    9   : ["R09", 0.2, 25, 8, [0.1340, 0.0190, 0.0070, 0.0000, 0.0000], [2500.0, 109.0, 30.0,    0.0,    0.0]],
    10  : ["R10", 1.0, 25, 8, [0.1380, 0.0920, 0.2300, 0.0920, 0.0000], [2000.0, 109.0, 30.0, 1500.0,    0.0]],
    # Logging Slash (light-logging-slash,medium-logging-slash,heavy-logging-slash)
    11  : ["R11", 1.0, 15, 8, [0.0690, 0.2070, 0.2530, 0.0000, 0.0000], [1500.0, 109.0, 30.0,    0.0,    0.0]],
    12  : ["R12", 2.3, 20, 8, [0.1840, 0.6440, 0.7590, 0.0000, 0.0000], [1500.0, 109.0, 30.0,    0.0,    0.0]],
    13  : ["R13", 3.0, 25, 8, [0.3220, 1.0580, 1.2880, 0.0000, 0.0000], [1500.0, 109.0, 30.0,    0.0,    0.0]],
    # Nonburnable (NB)
    91  : ["NB1", 0.0,  0, 0, [0.0000, 0.0000, 0.0000, 0.0000, 0.0000], [   0.0,   0.0,  0.0,    0.0,    0.0]],
    92  : ["NB2", 0.0,  0, 0, [0.0000, 0.0000, 0.0000, 0.0000, 0.0000], [   0.0,   0.0,  0.0,    0.0,    0.0]],
    93  : ["NB3", 0.0,  0, 0, [0.0000, 0.0000, 0.0000, 0.0000, 0.0000], [   0.0,   0.0,  0.0,    0.0,    0.0]],
    98  : ["NB4", 0.0,  0, 0, [0.0000, 0.0000, 0.0000, 0.0000, 0.0000], [   0.0,   0.0,  0.0,    0.0,    0.0]],
    99  : ["NB5", 0.0,  0, 0, [0.0000, 0.0000, 0.0000, 0.0000, 0.0000], [   0.0,   0.0,  0.0,    0.0,    0.0]],
    # Scott & Burgan 40:
    # Grass (GR)
    101 : ["GR1", 0.4, 15, 8, [0.0046, 0.0000, 0.0000, 0.0138, 0.0000], [2200.0, 109.0, 30.0, 2000.0,    0.0]],
    102 : ["GR2", 1.0, 15, 8, [0.0046, 0.0000, 0.0000, 0.0459, 0.0000], [2000.0, 109.0, 30.0, 1800.0,    0.0]],
    103 : ["GR3", 2.0, 30, 8, [0.0046, 0.0184, 0.0000, 0.0689, 0.0000], [1500.0, 109.0, 30.0, 1300.0,    0.0]],
    104 : ["GR4", 2.0, 15, 8, [0.0115, 0.0000, 0.0000, 0.0872, 0.0000], [2000.0, 109.0, 30.0, 1800.0,    0.0]],
    105 : ["GR5", 1.5, 40, 8, [0.0184, 0.0000, 0.0000, 0.1148, 0.0000], [1800.0, 109.0, 30.0, 1600.0,    0.0]],
    106 : ["GR6", 1.5, 40, 9, [0.0046, 0.0000, 0.0000, 0.1561, 0.0000], [2200.0, 109.0, 30.0, 2000.0,    0.0]],
    107 : ["GR7", 3.0, 15, 8, [0.0459, 0.0000, 0.0000, 0.2479, 0.0000], [2000.0, 109.0, 30.0, 1800.0,    0.0]],
    108 : ["GR8", 4.0, 30, 8, [0.0230, 0.0459, 0.0000, 0.3352, 0.0000], [1500.0, 109.0, 30.0, 1300.0,    0.0]],
    109 : ["GR9", 5.0, 40, 8, [0.0459, 0.0459, 0.0000, 0.4132, 0.0000], [1800.0, 109.0, 30.0, 1600.0,    0.0]],
    # Grass-Shrub (GS)
    121 : ["GS1", 0.9, 15, 8, [0.0092, 0.0000, 0.0000, 0.0230, 0.0298], [2000.0, 109.0, 30.0, 1800.0, 1800.0]],
    122 : ["GS2", 1.5, 15, 8, [0.0230, 0.0230, 0.0000, 0.0275, 0.0459], [2000.0, 109.0, 30.0, 1800.0, 1800.0]],
    123 : ["GS3", 1.8, 40, 8, [0.0138, 0.0115, 0.0000, 0.0666, 0.0574], [1800.0, 109.0, 30.0, 1600.0, 1600.0]],
    124 : ["GS4", 2.1, 40, 8, [0.0872, 0.0138, 0.0046, 0.1561, 0.3260], [1800.0, 109.0, 30.0, 1600.0, 1600.0]],
    # Shrub (SH)
    141 : ["SH1", 1.0, 15, 8, [0.0115, 0.0115, 0.0000, 0.0069, 0.0597], [2000.0, 109.0, 30.0, 1800.0, 1600.0]],
    142 : ["SH2", 1.0, 15, 8, [0.0620, 0.1102, 0.0344, 0.0000, 0.1768], [2000.0, 109.0, 30.0,    0.0, 1600.0]],
    143 : ["SH3", 2.4, 40, 8, [0.0207, 0.1377, 0.0000, 0.0000, 0.2847], [1600.0, 109.0, 30.0,    0.0, 1400.0]],
    144 : ["SH4", 3.0, 30, 8, [0.0390, 0.0528, 0.0092, 0.0000, 0.1171], [2000.0, 109.0, 30.0, 1800.0, 1600.0]],
    145 : ["SH5", 6.0, 15, 8, [0.1653, 0.0964, 0.0000, 0.0000, 0.1331], [ 750.0, 109.0, 30.0,    0.0, 1600.0]],
    146 : ["SH6", 2.0, 30, 8, [0.1331, 0.0666, 0.0000, 0.0000, 0.0643], [ 750.0, 109.0, 30.0,    0.0, 1600.0]],
    147 : ["SH7", 6.0, 15, 8, [0.1607, 0.2433, 0.1010, 0.0000, 0.1561], [ 750.0, 109.0, 30.0,    0.0, 1600.0]],
    148 : ["SH8", 3.0, 40, 8, [0.0941, 0.1561, 0.0390, 0.0000, 0.1997], [ 750.0, 109.0, 30.0,    0.0, 1600.0]],
    149 : ["SH9", 4.4, 40, 8, [0.2066, 0.1125, 0.0000, 0.0712, 0.3214], [ 750.0, 109.0, 30.0, 1800.0, 1500.0]],
    # Timber-Understory (TU)
    161 : ["TU1", 0.6, 20, 8, [0.0092, 0.0413, 0.0689, 0.0092, 0.0413], [2000.0, 109.0, 30.0, 1800.0, 1600.0]],
    162 : ["TU2", 1.0, 30, 8, [0.0436, 0.0826, 0.0574, 0.0000, 0.0092], [2000.0, 109.0, 30.0,    0.0, 1600.0]],
    163 : ["TU3", 1.3, 30, 8, [0.0505, 0.0069, 0.0115, 0.0298, 0.0505], [1800.0, 109.0, 30.0, 1600.0, 1400.0]],
    164 : ["TU4", 0.5, 12, 8, [0.2066, 0.0000, 0.0000, 0.0000, 0.0918], [2300.0, 109.0, 30.0,    0.0, 2000.0]],
    165 : ["TU5", 1.0, 25, 8, [0.1837, 0.1837, 0.1377, 0.0000, 0.1377], [1500.0, 109.0, 30.0,    0.0,  750.0]],
    # Timber Litter (TL)
    181 : ["TL1", 0.2, 30, 8, [0.0459, 0.1010, 0.1653, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    182 : ["TL2", 0.2, 25, 8, [0.0643, 0.1056, 0.1010, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    183 : ["TL3", 0.3, 20, 8, [0.0230, 0.1010, 0.1286, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    184 : ["TL4", 0.4, 25, 8, [0.0230, 0.0689, 0.1928, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    185 : ["TL5", 0.6, 25, 8, [0.0528, 0.1148, 0.2020, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0, 1600.0]],
    186 : ["TL6", 0.3, 25, 8, [0.1102, 0.0551, 0.0551, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    187 : ["TL7", 0.4, 25, 8, [0.0138, 0.0643, 0.3719, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    188 : ["TL8", 0.3, 35, 8, [0.2663, 0.0643, 0.0505, 0.0000, 0.0000], [1800.0, 109.0, 30.0,    0.0,    0.0]],
    189 : ["TL9", 0.6, 35, 8, [0.3053, 0.1515, 0.1905, 0.0000, 0.0000], [1800.0, 109.0, 30.0,    0.0, 1600.0]],
    # Slash-Blowdown (SB)
    201 : ["SB1", 1.0, 25, 8, [0.0689, 0.1377, 0.5051, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    202 : ["SB2", 1.0, 25, 8, [0.2066, 0.1951, 0.1837, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    203 : ["SB3", 1.2, 25, 8, [0.2525, 0.1263, 0.1377, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
    204 : ["SB4", 2.7, 25, 8, [0.2410, 0.1607, 0.2410, 0.0000, 0.0000], [2000.0, 109.0, 30.0,    0.0,    0.0]],
}
# fuel-model-compact-table ends here
# [[file:../../org/pyretechnics.org::expand-compact-fuel-model-table][expand-compact-fuel-model-table]]
def expand_compact_fuel_model(fuel_model_number):
    [name, delta, M_x_dead, h, w_o, sigma] = fuel_model_compact_table[fuel_model_number]
    [w_o_dead_1hr, w_o_dead_10hr, w_o_dead_100hr, w_o_live_herbaceous, w_o_live_woody] = w_o
    [sigma_dead_1hr, sigma_dead_10hr, sigma_dead_100hr, sigma_live_herbaceous, sigma_live_woody] = sigma
    M_x_dead = M_x_dead * 0.01
    h        = h * 1000.0
    # Conditionally set dead_herbaceous values
    dynamic               = fuel_model_number > 100 and w_o_live_herbaceous > 0.0
    M_x_dead_herbaceous   = M_x_dead              if dynamic else 0.0
    sigma_dead_herbaceous = sigma_live_herbaceous if dynamic else 0.0
    return {
        "name"    : name,
        "number"  : fuel_model_number,
        "delta"   : delta,
        "M_x"     : [M_x_dead, M_x_dead, M_x_dead, M_x_dead_herbaceous, 0.0, 0.0],
        "w_o"     : [w_o_dead_1hr, w_o_dead_10hr, w_o_dead_100hr, 0.0, w_o_live_herbaceous, w_o_live_woody],
        "sigma"   : [sigma_dead_1hr, sigma_dead_10hr, sigma_dead_100hr, sigma_dead_herbaceous, sigma_live_herbaceous, sigma_live_woody],
        "h"       : 6 * [h],
        "rho_p"   : 6 * [32.0],
        "S_T"     : 6 * [0.0555],
        "S_e"     : 6 * [0.01],
        "dynamic" : dynamic,
        "burnable": not (91 <= fuel_model_number <= 99),
    }


fuel_model_table = {k: expand_compact_fuel_model(k) for k in fuel_model_compact_table.keys()}

# expand-compact-fuel-model-table ends here
# [[file:../../org/pyretechnics.org::fuel-category-and-size-class-functions][fuel-category-and-size-class-functions]]
def map_category(f):
    return [f(0), f(1)]


def map_size_class(f):
    return [f(0), f(1), f(2), f(3), f(4), f(5)]


def category_sum(f):
    return f(0) + f(1)


def size_class_sum(f):
    return [f(0) + f(1) + f(2) + f(3), f(4) + f(5)]
# fuel-category-and-size-class-functions ends here
# [[file:../../org/pyretechnics.org::add-dynamic-fuel-loading][add-dynamic-fuel-loading]]
def add_dynamic_fuel_loading(fuel_model, M_f):
    if fuel_model["dynamic"]:
        # dynamic fuel model
        w_o                       = fuel_model["w_o"]
        live_herbaceous_load      = w_o[4]
        live_herbaceous_moisture  = M_f[4]
        fraction_green            = max(0.0, min(1.0, (live_herbaceous_moisture / 0.9) - 0.3333333333333333))
        fraction_cured            = 1.0 - fraction_green
        dynamic_fuel_model        = fuel_model.copy() # shallow copy
        dynamic_fuel_model["M_f"] = [
            M_f[0],
            M_f[1],
            M_f[2],
            M_f[0], # set dead_herbaceous to dead_1hr
            M_f[4],
            M_f[5],
        ]
        dynamic_fuel_model["w_o"] = [
            w_o[0],
            w_o[1],
            w_o[2],
            live_herbaceous_load * fraction_cured, # dead_herbaceous
            live_herbaceous_load * fraction_green, # live_herbaceous
            w_o[5],
        ]
        return dynamic_fuel_model
    else:
        # static fuel model
        static_fuel_model = fuel_model.copy() # shallow copy
        static_fuel_model["M_f"] = M_f
        return static_fuel_model
# add-dynamic-fuel-loading ends here
# [[file:../../org/pyretechnics.org::add-weighting-factors][add-weighting-factors]]
def add_weighting_factors(fuel_model):
    w_o                         = fuel_model["w_o"]
    sigma                       = fuel_model["sigma"]
    rho_p                       = fuel_model["rho_p"]
    def msc_Aij(i):
        return (sigma[i] * w_o[i]) / rho_p[i]
    A_ij                        = map_size_class(msc_Aij)
    def scs_A_ij(i):
        return A_ij[i]
    A_i                         = size_class_sum(scs_A_ij)
    def scs_A_i(i):
        return A_i[i]
    A_T                         = category_sum(scs_A_i)
    def msc_fij(i):
        A = A_i[i//4] # FIXME clever but unclear
        return (A_ij[i] / A) if (A > 0.0) else 0.0
    f_ij                        = map_size_class(msc_fij)
    def msc_f_i(i):
        return (A_i[i] / A_T) if A_T > 0.0 else 0.0
    f_i                         = map_category(msc_f_i)
    def msc_firemod_size_class(i):
        s = sigma[i]
        return (
            1 if (s >= 1200.0)
            else 2 if (s >= 192.0)
            else 3 if (s >= 96.0)
            else 4 if (s >= 48.0)
            else 5 if (s >= 16.0)
            else 6
        )
    firemod_size_classes        = map_size_class(msc_firemod_size_class)
    def msc_gij(i):
        c = firemod_size_classes[i]
        return ((
            (f_ij[0] if (c == firemod_size_classes[0]) else 0.0)
            + (f_ij[1] if (c == firemod_size_classes[1]) else 0.0)
                + (f_ij[2] if (c == firemod_size_classes[2]) else 0.0)
                + (f_ij[3] if (c == firemod_size_classes[3]) else 0.0))
            if (i < 4) else
            ((f_ij[4] if (c == firemod_size_classes[4]) else 0.0)
                + (f_ij[5] if (c == firemod_size_classes[5]) else 0.0)))
            
    g_ij                        = map_size_class(msc_gij)
    weighted_fuel_model         = fuel_model.copy() # shallow copy
    weighted_fuel_model["f_ij"] = f_ij
    weighted_fuel_model["f_i"]  = f_i
    weighted_fuel_model["g_ij"] = g_ij
    return weighted_fuel_model
# add-weighting-factors ends here
# [[file:../../org/pyretechnics.org::add-live-moisture-of-extinction][add-live-moisture-of-extinction]]
from math import exp

def add_live_moisture_of_extinction(fuel_model):
    """
    Equation 88 from Rothermel 1972 adjusted by Albini 1976 Appendix III.
    """
    w_o                       = fuel_model["w_o"]
    sigma                     = fuel_model["sigma"]
    M_f                       = fuel_model["M_f"]
    M_x                       = fuel_model["M_x"]
    def msc_loading_factor(i):
        sigma_ij = sigma[i]
        A = -138.0 if (i < 4) else -500.0
        return w_o[i] * exp(A / sigma_ij) if (sigma_ij > 0.0) else 0.0
    loading_factors           = map_size_class(msc_loading_factor)
    def scs_loading_factor(i):
        return loading_factors[i]
    [dead_loading_factor,
     live_loading_factor]     = size_class_sum(scs_loading_factor)
    def scs_dead_moisture_factor(i):
        return M_f[i] * loading_factors[i]
    [dead_moisture_factor, _] = size_class_sum(scs_dead_moisture_factor)
    dead_to_live_ratio        = (dead_loading_factor / live_loading_factor) if (live_loading_factor > 0.0) else None
    dead_fuel_moisture        = (dead_moisture_factor / dead_loading_factor) if (dead_loading_factor > 0.0) else 0.0
    M_x_dead                  = M_x[0]
    M_x_live                  = max(M_x_dead,
                                    (2.9 * dead_to_live_ratio * (1.0 - (dead_fuel_moisture / M_x_dead))) - 0.226
                                    ) if (live_loading_factor > 0.0) else M_x_dead
    moisturized_fuel_model    = fuel_model.copy() # shallow copy
    moisturized_fuel_model["M_x"] = [
        M_x[0],
        M_x[1],
        M_x[2],
        M_x[3],
        M_x_live,
        M_x_live,
    ]
    return moisturized_fuel_model
# add-live-moisture-of-extinction ends here
# [[file:../../org/pyretechnics.org::moisturize][moisturize]]
# TODO: If these functions aren't called anywhere else, create a copy
#       of the fuel model here and mutate it in the called functions.
def moisturize(fuel_model, fuel_moisture):
    dynamic_fuel_model     = add_dynamic_fuel_loading(fuel_model, fuel_moisture)
    weighted_fuel_model    = add_weighting_factors(dynamic_fuel_model)
    moisturized_fuel_model = add_live_moisture_of_extinction(weighted_fuel_model)
    return moisturized_fuel_model
# moisturize ends here
