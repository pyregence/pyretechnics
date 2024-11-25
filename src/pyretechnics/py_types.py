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
# py-types-py ends here
