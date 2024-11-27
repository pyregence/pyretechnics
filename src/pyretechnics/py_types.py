# [[file:../../org/pyretechnics.org::py-types-py][py-types-py]]
import cython as cy
from 

#==============================================================
# Runtime-defined type aliases
#==============================================================

# TODO: Maybe pyidx should be cy.int?
# TODO: Maybe we should use C arrays instead of tuples(structs)?
pyidx     = cy.typedef(cy.Py_ssize_t)
vec_xy    = cy.typedef(tuple[cy.float, cy.float])
vec_xyz   = cy.typedef(tuple[cy.float, cy.float, cy.float])
coord_yx  = cy.typedef(tuple[pyidx, pyidx])
coord_tyx = cy.typedef(tuple[pyidx, pyidx, pyidx])
fcatarr   = cy.typedef(tuple[cy.float, cy.float])
fclaarr   = cy.typedef(tuple[cy.float, cy.float, cy.float, cy.float, cy.float, cy.float])

FireBehaviorMax = cy.struct(
    max_fire_type         : cy.int,
    max_spread_rate       : cy.float,
    max_spread_direction  : vec_xyz,
    max_fireline_intensity: cy.float,
    max_flame_length      : cy.float,
    length_to_width_ratio : cy.float,
    eccentricity          : cy.float,
    critical_spread_rate  : cy.float
)

SpreadBehavior = cy.struct(
    dphi_dt            : cy.float,
    fire_type          : cy.int,
    spread_rate        : cy.float,
    spread_direction   : vec_xyz,
    fireline_intensity : cy.float,
    flame_length       : cy.float,
)

FuelModel = cy.struct(
    # name    = cy.string,
    number   : cy.int,
    delta    : cy.float,
    M_x      : fclaarr,
    M_f      : fclaarr,
    w_o      : fclaarr,
    sigma    : fclaarr,
    h        : fclaarr,
    rho_p    : fclaarr,
    S_T      : fclaarr,
    S_e      : fclaarr,
    dynamic  : cy.bint,
    burnable : cy.bint,
    f_ij     : fclaarr
    f_i      : fcatarr
    g_ij     : fclaarr
)
# py-types-py ends here
