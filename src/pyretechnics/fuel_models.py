# [[file:../../org/pyretechnics.org::fuel-models-imports][fuel-models-imports]]
import cython
import cython as cy
if cython.compiled:
    from cython.cimports.libc.math import exp
    from cython.cimports.pyretechnics.cy_types import fcatarr, fclaarr, CompactFuelModel, FuelModel
else:
    from math import exp
    from pyretechnics.py_types import fcatarr, fclaarr, CompactFuelModel, FuelModel
# fuel-models-imports ends here
# [[file:../../org/pyretechnics.org::fuel-model-compact-table][fuel-model-compact-table]]
# TODO: OPTIM Replace this dictionary with something more efficient
# Lookup table including entries for each of the Anderson 13 and Scott & Burgan 40 fuel models.
#
# The fields have the following meanings:
#   {
#     fuel_model_number: (delta,
#                         M_x_dead,
#                         h,
#                         w_o_dead_1hr,
#                         w_o_dead_10hr,
#                         w_o_dead_100hr,
#                         w_o_live_herbaceous,
#                         w_o_live_woody,
#                         sigma_dead_1hr,
#                         sigma_dead_10hr,
#                         sigma_dead_100hr,
#                         sigma_live_herbaceous,
#                         sigma_live_woody), # name
#   }
compact_fuel_model_table = cy.declare(dict, { # dict[int, CompactFuelModel]
    # Anderson 13:
    # Grass and Grass-dominated (short-grass,timber-grass-and-understory,tall-grass)
    1  : (1.0, 12.0, 8.0, 0.0340, 0.0000, 0.0000, 0.0000, 0.0000, 3500.0,   0.0,  0.0,    0.0,    0.0), # R01
    2  : (1.0, 15.0, 8.0, 0.0920, 0.0460, 0.0230, 0.0230, 0.0000, 3000.0, 109.0, 30.0, 1500.0,    0.0), # R02
    3  : (2.5, 25.0, 8.0, 0.1380, 0.0000, 0.0000, 0.0000, 0.0000, 1500.0,   0.0,  0.0,    0.0,    0.0), # R03
    # Chaparral and Shrubfields (chaparral,brush,dormant-brush-hardwood-slash,southern-rough)
    4  : (6.0, 20.0, 8.0, 0.2300, 0.1840, 0.0920, 0.2300, 0.0000, 2000.0, 109.0, 30.0, 1500.0,    0.0), # R04
    5  : (2.0, 20.0, 8.0, 0.0460, 0.0230, 0.0000, 0.0920, 0.0000, 2000.0, 109.0,  0.0, 1500.0,    0.0), # R05
    6  : (2.5, 25.0, 8.0, 0.0690, 0.1150, 0.0920, 0.0000, 0.0000, 1750.0, 109.0, 30.0,    0.0,    0.0), # R06
    7  : (2.5, 40.0, 8.0, 0.0520, 0.0860, 0.0690, 0.0170, 0.0000, 1750.0, 109.0, 30.0, 1550.0,    0.0), # R07
    # Timber Litter (closed-timber-litter,hardwood-litter,timber-litter-and-understory)
    8  : (0.2, 30.0, 8.0, 0.0690, 0.0460, 0.1150, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # R08
    9  : (0.2, 25.0, 8.0, 0.1340, 0.0190, 0.0070, 0.0000, 0.0000, 2500.0, 109.0, 30.0,    0.0,    0.0), # R09
    10 : (1.0, 25.0, 8.0, 0.1380, 0.0920, 0.2300, 0.0920, 0.0000, 2000.0, 109.0, 30.0, 1500.0,    0.0), # R10
    # Logging Slash (light-logging-slash,medium-logging-slash,heavy-logging-slash)
    11 : (1.0, 15.0, 8.0, 0.0690, 0.2070, 0.2530, 0.0000, 0.0000, 1500.0, 109.0, 30.0,    0.0,    0.0), # R11
    12 : (2.3, 20.0, 8.0, 0.1840, 0.6440, 0.7590, 0.0000, 0.0000, 1500.0, 109.0, 30.0,    0.0,    0.0), # R12
    13 : (3.0, 25.0, 8.0, 0.3220, 1.0580, 1.2880, 0.0000, 0.0000, 1500.0, 109.0, 30.0,    0.0,    0.0), # R13
    # Nonburnable (NB)
    91 : (0.0,  0.0, 0.0, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,    0.0,   0.0,  0.0,    0.0,    0.0), # NB1
    92 : (0.0,  0.0, 0.0, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,    0.0,   0.0,  0.0,    0.0,    0.0), # NB2
    93 : (0.0,  0.0, 0.0, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,    0.0,   0.0,  0.0,    0.0,    0.0), # NB3
    98 : (0.0,  0.0, 0.0, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,    0.0,   0.0,  0.0,    0.0,    0.0), # NB4
    99 : (0.0,  0.0, 0.0, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,    0.0,   0.0,  0.0,    0.0,    0.0), # NB5
    # Scott & Burgan 40:
    # Grass (GR)
    101: (0.4, 15.0, 8.0, 0.0046, 0.0000, 0.0000, 0.0138, 0.0000, 2200.0, 109.0, 30.0, 2000.0,    0.0), # GR1
    102: (1.0, 15.0, 8.0, 0.0046, 0.0000, 0.0000, 0.0459, 0.0000, 2000.0, 109.0, 30.0, 1800.0,    0.0), # GR2
    103: (2.0, 30.0, 8.0, 0.0046, 0.0184, 0.0000, 0.0689, 0.0000, 1500.0, 109.0, 30.0, 1300.0,    0.0), # GR3
    104: (2.0, 15.0, 8.0, 0.0115, 0.0000, 0.0000, 0.0872, 0.0000, 2000.0, 109.0, 30.0, 1800.0,    0.0), # GR4
    105: (1.5, 40.0, 8.0, 0.0184, 0.0000, 0.0000, 0.1148, 0.0000, 1800.0, 109.0, 30.0, 1600.0,    0.0), # GR5
    106: (1.5, 40.0, 9.0, 0.0046, 0.0000, 0.0000, 0.1561, 0.0000, 2200.0, 109.0, 30.0, 2000.0,    0.0), # GR6
    107: (3.0, 15.0, 8.0, 0.0459, 0.0000, 0.0000, 0.2479, 0.0000, 2000.0, 109.0, 30.0, 1800.0,    0.0), # GR7
    108: (4.0, 30.0, 8.0, 0.0230, 0.0459, 0.0000, 0.3352, 0.0000, 1500.0, 109.0, 30.0, 1300.0,    0.0), # GR8
    109: (5.0, 40.0, 8.0, 0.0459, 0.0459, 0.0000, 0.4132, 0.0000, 1800.0, 109.0, 30.0, 1600.0,    0.0), # GR9
    # Grass-Shrub (GS)
    121: (0.9, 15.0, 8.0, 0.0092, 0.0000, 0.0000, 0.0230, 0.0298, 2000.0, 109.0, 30.0, 1800.0, 1800.0), # GS1
    122: (1.5, 15.0, 8.0, 0.0230, 0.0230, 0.0000, 0.0275, 0.0459, 2000.0, 109.0, 30.0, 1800.0, 1800.0), # GS2
    123: (1.8, 40.0, 8.0, 0.0138, 0.0115, 0.0000, 0.0666, 0.0574, 1800.0, 109.0, 30.0, 1600.0, 1600.0), # GS3
    124: (2.1, 40.0, 8.0, 0.0872, 0.0138, 0.0046, 0.1561, 0.3260, 1800.0, 109.0, 30.0, 1600.0, 1600.0), # GS4
    # Shrub (SH)
    141: (1.0, 15.0, 8.0, 0.0115, 0.0115, 0.0000, 0.0069, 0.0597, 2000.0, 109.0, 30.0, 1800.0, 1600.0), # SH1
    142: (1.0, 15.0, 8.0, 0.0620, 0.1102, 0.0344, 0.0000, 0.1768, 2000.0, 109.0, 30.0,    0.0, 1600.0), # SH2
    143: (2.4, 40.0, 8.0, 0.0207, 0.1377, 0.0000, 0.0000, 0.2847, 1600.0, 109.0, 30.0,    0.0, 1400.0), # SH3
    144: (3.0, 30.0, 8.0, 0.0390, 0.0528, 0.0092, 0.0000, 0.1171, 2000.0, 109.0, 30.0, 1800.0, 1600.0), # SH4
    145: (6.0, 15.0, 8.0, 0.1653, 0.0964, 0.0000, 0.0000, 0.1331,  750.0, 109.0, 30.0,    0.0, 1600.0), # SH5
    146: (2.0, 30.0, 8.0, 0.1331, 0.0666, 0.0000, 0.0000, 0.0643,  750.0, 109.0, 30.0,    0.0, 1600.0), # SH6
    147: (6.0, 15.0, 8.0, 0.1607, 0.2433, 0.1010, 0.0000, 0.1561,  750.0, 109.0, 30.0,    0.0, 1600.0), # SH7
    148: (3.0, 40.0, 8.0, 0.0941, 0.1561, 0.0390, 0.0000, 0.1997,  750.0, 109.0, 30.0,    0.0, 1600.0), # SH8
    149: (4.4, 40.0, 8.0, 0.2066, 0.1125, 0.0000, 0.0712, 0.3214,  750.0, 109.0, 30.0, 1800.0, 1500.0), # SH9
    # Timber-Understory (TU)
    161: (0.6, 20.0, 8.0, 0.0092, 0.0413, 0.0689, 0.0092, 0.0413, 2000.0, 109.0, 30.0, 1800.0, 1600.0), # TU1
    162: (1.0, 30.0, 8.0, 0.0436, 0.0826, 0.0574, 0.0000, 0.0092, 2000.0, 109.0, 30.0,    0.0, 1600.0), # TU2
    163: (1.3, 30.0, 8.0, 0.0505, 0.0069, 0.0115, 0.0298, 0.0505, 1800.0, 109.0, 30.0, 1600.0, 1400.0), # TU3
    164: (0.5, 12.0, 8.0, 0.2066, 0.0000, 0.0000, 0.0000, 0.0918, 2300.0, 109.0, 30.0,    0.0, 2000.0), # TU4
    165: (1.0, 25.0, 8.0, 0.1837, 0.1837, 0.1377, 0.0000, 0.1377, 1500.0, 109.0, 30.0,    0.0,  750.0), # TU5
    # Timber Litter (TL)
    181: (0.2, 30.0, 8.0, 0.0459, 0.1010, 0.1653, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # TL1
    182: (0.2, 25.0, 8.0, 0.0643, 0.1056, 0.1010, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # TL2
    183: (0.3, 20.0, 8.0, 0.0230, 0.1010, 0.1286, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # TL3
    184: (0.4, 25.0, 8.0, 0.0230, 0.0689, 0.1928, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # TL4
    185: (0.6, 25.0, 8.0, 0.0528, 0.1148, 0.2020, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0, 1600.0), # TL5
    186: (0.3, 25.0, 8.0, 0.1102, 0.0551, 0.0551, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # TL6
    187: (0.4, 25.0, 8.0, 0.0138, 0.0643, 0.3719, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # TL7
    188: (0.3, 35.0, 8.0, 0.2663, 0.0643, 0.0505, 0.0000, 0.0000, 1800.0, 109.0, 30.0,    0.0,    0.0), # TL8
    189: (0.6, 35.0, 8.0, 0.3053, 0.1515, 0.1905, 0.0000, 0.0000, 1800.0, 109.0, 30.0,    0.0, 1600.0), # TL9
    # Slash-Blowdown (SB)
    201: (1.0, 25.0, 8.0, 0.0689, 0.1377, 0.5051, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # SB1
    202: (1.0, 25.0, 8.0, 0.2066, 0.1951, 0.1837, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # SB2
    203: (1.2, 25.0, 8.0, 0.2525, 0.1263, 0.1377, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # SB3
    204: (2.7, 25.0, 8.0, 0.2410, 0.1607, 0.2410, 0.0000, 0.0000, 2000.0, 109.0, 30.0,    0.0,    0.0), # SB4
})
# fuel-model-compact-table ends here
# [[file:../../org/pyretechnics.org::expand-compact-fuel-model-table][expand-compact-fuel-model-table]]
# NOTE: We use this variable for X > 0.0 comparisons to account for floating point precision issues.
almost_zero = cy.declare(cy.float, 1e-6)


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def is_burnable_fuel_model_number(fuel_model_number: cy.int) -> cy.bint:
    return not (91 <= fuel_model_number <= 99)


@cy.cfunc
@cy.exceptval(check=False)
def compute_exp_A_sigma(A: cy.float, sigma_ij: cy.float) -> cy.float:
    if sigma_ij > almost_zero:
        return exp(A / sigma_ij)
    else:
        return 0.0


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def compute_firemod_size_class(sigma_i: cy.float) -> cy.float:
    return (
        1.0 if (sigma_i >= 1200.0)
        else 2.0 if (sigma_i >= 192.0)
        else 3.0 if (sigma_i >= 96.0)
        else 4.0 if (sigma_i >= 48.0)
        else 5.0 if (sigma_i >= 16.0)
        else 6.0
    )


@cy.cfunc
@cy.exceptval(check=False)
def expand_compact_fuel_model(fuel_model_number: cy.int) -> FuelModel:
    # Look up the CompactFuelModel by fuel_model_number
    cfm: CompactFuelModel = compact_fuel_model_table[fuel_model_number]
    # Unpack the CompactFuelModel values
    delta                : cy.float = cfm[0]
    M_x_dead             : cy.float = cfm[1]
    h                    : cy.float = cfm[2]
    w_o_dead_1hr         : cy.float = cfm[3]
    w_o_dead_10hr        : cy.float = cfm[4]
    w_o_dead_100hr       : cy.float = cfm[5]
    w_o_live_herbaceous  : cy.float = cfm[6]
    w_o_live_woody       : cy.float = cfm[7]
    sigma_dead_1hr       : cy.float = cfm[8]
    sigma_dead_10hr      : cy.float = cfm[9]
    sigma_dead_100hr     : cy.float = cfm[10]
    sigma_live_herbaceous: cy.float = cfm[11]
    sigma_live_woody     : cy.float = cfm[12]
    # Expand compressed values
    M_x_dead: cy.float = M_x_dead * 0.01
    h       : cy.float = h * 1000.0
    # Pre-compute some dynamic fuel model values
    dynamic              : cy.bint  = fuel_model_number > 100 and w_o_live_herbaceous > almost_zero
    M_x_dead_herbaceous  : cy.float = M_x_dead              if dynamic else 0.0
    sigma_dead_herbaceous: cy.float = sigma_live_herbaceous if dynamic else 0.0
    # Re-pack everything into a FuelModel struct
    return FuelModel(
        number               = fuel_model_number,
        delta                = delta,
        M_x                  = (M_x_dead,
                                M_x_dead,
                                M_x_dead,
                                M_x_dead_herbaceous,
                                0.0,
                                0.0),
        M_f                  = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        w_o                  = (w_o_dead_1hr,
                                w_o_dead_10hr,
                                w_o_dead_100hr,
                                0.0,
                                w_o_live_herbaceous,
                                w_o_live_woody),
        sigma                = (sigma_dead_1hr,
                                sigma_dead_10hr,
                                sigma_dead_100hr,
                                sigma_dead_herbaceous,
                                sigma_live_herbaceous,
                                sigma_live_woody),
        h                    = (h, h, h, h, h, h),
        rho_p                = (32.0, 32.0, 32.0, 32.0, 32.0, 32.0),
        S_T                  = (0.0555, 0.0555, 0.0555, 0.0555, 0.0555, 0.0555),
        S_e                  = (0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
        dynamic              = dynamic,
        burnable             = is_burnable_fuel_model_number(fuel_model_number),
        exp_A_sigma          = (compute_exp_A_sigma(-138.0, sigma_dead_1hr),
                                compute_exp_A_sigma(-138.0, sigma_dead_10hr),
                                compute_exp_A_sigma(-138.0, sigma_dead_100hr),
                                compute_exp_A_sigma(-138.0, sigma_dead_herbaceous),
                                compute_exp_A_sigma(-500.0, sigma_live_herbaceous),
                                compute_exp_A_sigma(-500.0, sigma_live_woody)),
        firemod_size_classes = (compute_firemod_size_class(sigma_dead_1hr),
                                compute_firemod_size_class(sigma_dead_10hr),
                                compute_firemod_size_class(sigma_dead_100hr),
                                compute_firemod_size_class(sigma_dead_herbaceous),
                                compute_firemod_size_class(sigma_live_herbaceous),
                                compute_firemod_size_class(sigma_live_woody)),
        f_ij                 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        f_i                  = (0.0, 0.0),
        g_ij                 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    )


# TODO: OPTIM Replace this dictionary with something more efficient
fuel_model_table = cy.declare(dict, { # dict[int, FuelModel]
    k: expand_compact_fuel_model(k) for k in compact_fuel_model_table.keys()
})


@cy.ccall
@cy.inline
def list_fuel_model_numbers() -> list[int]:
    return list(fuel_model_table.keys())


@cy.ccall
@cy.inline
def list_fuel_models() -> list[FuelModel]:
    return list(fuel_model_table.values())


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def fuel_model_exists(fuel_model_number: cy.int) -> cy.bint:
    return fuel_model_number in fuel_model_table


@cy.ccall
@cy.inline
@cy.exceptval(check=False)
def get_fuel_model(fuel_model_number: cy.int) -> FuelModel:
    return cy.cast(FuelModel, fuel_model_table[fuel_model_number])
# expand-compact-fuel-model-table ends here
# [[file:../../org/pyretechnics.org::add-dynamic-fuel-loading][add-dynamic-fuel-loading]]
@cy.cfunc
@cy.exceptval(check=False)
def add_dynamic_fuel_loading(fuel_model: FuelModel, M_f: fclaarr) -> FuelModel:
    """
    Updates M_f and w_o. Mutates and returns fuel_model.
    """
    if fuel_model.dynamic:
        # === Dynamic Fuel Model ===
        # Calculate fraction_green and fraction_cured
        w_o                     : fclaarr  = fuel_model.w_o
        live_herbaceous_load    : cy.float = w_o[4]
        live_herbaceous_moisture: cy.float = M_f[4]
        fraction_green          : cy.float = max(0.0, min(1.0, (live_herbaceous_moisture / 0.9) - 0.3333333333333333))
        fraction_cured          : cy.float = 1.0 - fraction_green
        # Set M_f[dead_herbaceous] equal to M_f[dead_1hr]
        M_f_dynamic: fclaarr = (M_f[0], # dead_1hr
                                M_f[1],
                                M_f[2],
                                M_f[0], # dead_herbaceous
                                M_f[4],
                                M_f[5])
        # Set w_o[dead_herbaceous] and w_o[live_herbaceous] based on fraction_cured and fraction_green
        w_o_dynamic: fclaarr = (w_o[0],
                                w_o[1],
                                w_o[2],
                                live_herbaceous_load * fraction_cured, # dead_herbaceous
                                live_herbaceous_load * fraction_green, # live_herbaceous
                                w_o[5])
        # Update the passed in fuel_model to use M_f_dynamic and w_o_dynamic
        fuel_model.M_f = M_f_dynamic
        fuel_model.w_o = w_o_dynamic
        return fuel_model
    else:
        # === Static Fuel Model ===
        # Update the passed in fuel_model to use the passed in M_f
        fuel_model.M_f = M_f
        return fuel_model
# add-dynamic-fuel-loading ends here
# [[file:../../org/pyretechnics.org::add-weighting-factors][add-weighting-factors]]
# TODO: OPTIM pre-compute this conditional branching since it's fully determined by sigma.
#       This information might be represented efficiently in bit flags.
@cy.cfunc
@cy.exceptval(check=False)
def compute_gij(firemod_size_classes : fclaarr,
                f_ij                 : fclaarr,
                firemod_size_class_ij: cy.float,
                is_dead              : cy.bint) -> cy.float:
    """
    Sums the f_ij of the same category (dead/live) as i, and having the same firemod_size_class.

    NOTE: There may be repetitions in firemod_size_classes, which is why this expression is not
          trivially equal to f_ij[i].
    """
    if is_dead:
        f_ij_0: cy.float = (f_ij[0] if (firemod_size_class_ij == firemod_size_classes[0]) else 0.0)
        f_ij_1: cy.float = (f_ij[1] if (firemod_size_class_ij == firemod_size_classes[1]) else 0.0)
        f_ij_2: cy.float = (f_ij[2] if (firemod_size_class_ij == firemod_size_classes[2]) else 0.0)
        f_ij_3: cy.float = (f_ij[3] if (firemod_size_class_ij == firemod_size_classes[3]) else 0.0)
        return f_ij_0 + f_ij_1 + f_ij_2 + f_ij_3
    else:
        f_ij_4: cy.float = (f_ij[4] if (firemod_size_class_ij == firemod_size_classes[4]) else 0.0)
        f_ij_5: cy.float = (f_ij[5] if (firemod_size_class_ij == firemod_size_classes[5]) else 0.0)
        return f_ij_4 + f_ij_5


@cy.cfunc
@cy.exceptval(check=False)
def add_weighting_factors(fuel_model: FuelModel) -> FuelModel:
    """
    Assigns f_ij, f_i, and g_ij. Mutates and returns fuel_model.
    """
    # Calculate A_ij, A_i, and A_T
    w_o  : fclaarr  = fuel_model.w_o
    sigma: fclaarr  = fuel_model.sigma
    rho_p: fclaarr  = fuel_model.rho_p
    A_ij : fclaarr  = ((sigma[0] * w_o[0]) / rho_p[0],
                       (sigma[1] * w_o[1]) / rho_p[1],
                       (sigma[2] * w_o[2]) / rho_p[2],
                       (sigma[3] * w_o[3]) / rho_p[3],
                       (sigma[4] * w_o[4]) / rho_p[4],
                       (sigma[5] * w_o[5]) / rho_p[5]) # TODO: OPTIM pre-compute sigma/rho_p
    A_i_0: cy.float = A_ij[0] + A_ij[1] + A_ij[2] + A_ij[3]
    A_i_1: cy.float = A_ij[4] + A_ij[5]
    A_i  : fcatarr  = (A_i_0, A_i_1)
    A_T  : cy.float = A_i_0 + A_i_1
    # Calculate f_ij
    f_ij_0   : cy.float = 0.0
    f_ij_1   : cy.float = 0.0
    f_ij_2   : cy.float = 0.0
    f_ij_3   : cy.float = 0.0
    f_ij_4   : cy.float = 0.0
    f_ij_5   : cy.float = 0.0
    A_i_0_inv: cy.float
    A_i_1_inv: cy.float
    if A_i_0 > almost_zero:
        A_i_0_inv = 1.0 / A_i_0
        f_ij_0    = A_ij[0] * A_i_0_inv
        f_ij_1    = A_ij[1] * A_i_0_inv
        f_ij_2    = A_ij[2] * A_i_0_inv
        f_ij_3    = A_ij[3] * A_i_0_inv
    if A_i_1 > almost_zero:
        A_i_1_inv = 1.0 / A_i_1
        f_ij_4    = A_ij[4] * A_i_1_inv
        f_ij_5    = A_ij[5] * A_i_1_inv
    f_ij: fclaarr = (f_ij_0, f_ij_1, f_ij_2, f_ij_3, f_ij_4, f_ij_5)
    # Calculate f_i
    f_i_0  : cy.float = 0.0
    f_i_1  : cy.float = 0.0
    A_T_inv: cy.float
    if A_T > almost_zero:
        A_T_inv = 1.0 / A_T
        f_i_0   = A_i_0 * A_T_inv
        f_i_1   = A_i_1 * A_T_inv
    f_i: fcatarr = (f_i_0, f_i_1)
    # Calculate g_ij
    firemod_size_classes: fclaarr = fuel_model.firemod_size_classes
    g_ij                : fclaarr = (compute_gij(firemod_size_classes, f_ij, firemod_size_classes[0], True),
                                     compute_gij(firemod_size_classes, f_ij, firemod_size_classes[1], True),
                                     compute_gij(firemod_size_classes, f_ij, firemod_size_classes[2], True),
                                     compute_gij(firemod_size_classes, f_ij, firemod_size_classes[3], True),
                                     compute_gij(firemod_size_classes, f_ij, firemod_size_classes[4], False),
                                     compute_gij(firemod_size_classes, f_ij, firemod_size_classes[5], False))
    # Update the passed in fuel_model to use f_ij, f_i, and g_ij
    fuel_model.f_ij = f_ij
    fuel_model.f_i  = f_i
    fuel_model.g_ij = g_ij
    return fuel_model
# add-weighting-factors ends here
# [[file:../../org/pyretechnics.org::add-live-moisture-of-extinction][add-live-moisture-of-extinction]]
@cy.cfunc
@cy.exceptval(check=False)
def add_live_moisture_of_extinction(fuel_model: FuelModel) -> FuelModel:
    """
    Equation 88 from Rothermel 1972 adjusted by Albini 1976 Appendix III.

    Updates M_x. Mutates and returns fuel_model.
    """
    # Calculate dead_moisture_factor, dead_loading_factor, and live_loading_factor
    w_o                 : fclaarr  = fuel_model.w_o
    exp_A_sigma         : fclaarr  = fuel_model.exp_A_sigma
    M_f                 : fclaarr  = fuel_model.M_f
    M_x                 : fclaarr  = fuel_model.M_x
    M_x_dead            : cy.float = M_x[0]
    loading_factors     : fclaarr  = (w_o[0] * exp_A_sigma[0],
                                      w_o[1] * exp_A_sigma[1],
                                      w_o[2] * exp_A_sigma[2],
                                      w_o[3] * exp_A_sigma[3],
                                      w_o[4] * exp_A_sigma[4],
                                      w_o[5] * exp_A_sigma[5])
    dead_moisture_factor: cy.float = (M_f[0] * loading_factors[0] +
                                      M_f[1] * loading_factors[1] +
                                      M_f[2] * loading_factors[2] +
                                      M_f[3] * loading_factors[3])
    dead_loading_factor : cy.float = loading_factors[0] + loading_factors[1] + loading_factors[2] + loading_factors[3]
    live_loading_factor : cy.float = loading_factors[4] + loading_factors[5]
    # Calculate M_x_live
    dead_fuel_moisture: cy.float
    dead_to_live_ratio: cy.float
    M_x_live          : cy.float
    if (dead_loading_factor > almost_zero and live_loading_factor > almost_zero):
        dead_fuel_moisture = dead_moisture_factor / dead_loading_factor
        dead_to_live_ratio = dead_loading_factor / live_loading_factor
        M_x_live           = max(M_x_dead,
                                 (2.9 * dead_to_live_ratio * (1.0 - (dead_fuel_moisture / M_x_dead))) - 0.226)
    else:
        M_x_live = M_x_dead
    # Calculate M_x_new
    M_x_new: fclaarr = (M_x[0],
                        M_x[1],
                        M_x[2],
                        M_x[3],
                        M_x_live,
                        M_x_live)
    # Update the passed in fuel_model to use M_x_new
    fuel_model.M_x = M_x_new
    return fuel_model
# add-live-moisture-of-extinction ends here
# [[file:../../org/pyretechnics.org::moisturize][moisturize]]
@cy.ccall
@cy.exceptval(check=False)
def moisturize(fuel_model: FuelModel, fuel_moisture: fclaarr) -> FuelModel:
    """
    Updates w_o, M_f, and M_x and assigns f_ij, f_i, and g_ij.
    Returns a new FuelModel struct.
    """
    # Create a copy of fuel_model
    fuel_model_copy: FuelModel = FuelModel(
        number               = fuel_model.number,
        delta                = fuel_model.delta,
        M_x                  = fuel_model.M_x,
        M_f                  = fuel_model.M_f,
        w_o                  = fuel_model.w_o,
        sigma                = fuel_model.sigma,
        h                    = fuel_model.h,
        rho_p                = fuel_model.rho_p,
        S_T                  = fuel_model.S_T,
        S_e                  = fuel_model.S_e,
        dynamic              = fuel_model.dynamic,
        burnable             = fuel_model.burnable,
        exp_A_sigma          = fuel_model.exp_A_sigma,
        firemod_size_classes = fuel_model.firemod_size_classes,
        f_ij                 = fuel_model.f_ij,
        f_i                  = fuel_model.f_i,
        g_ij                 = fuel_model.g_ij,
    )
    # Mutate and return the copy
    dynamic_fuel_model    : FuelModel = add_dynamic_fuel_loading(fuel_model_copy, fuel_moisture)
    weighted_fuel_model   : FuelModel = add_weighting_factors(dynamic_fuel_model)
    moisturized_fuel_model: FuelModel = add_live_moisture_of_extinction(weighted_fuel_model)
    return moisturized_fuel_model
# moisturize ends here
