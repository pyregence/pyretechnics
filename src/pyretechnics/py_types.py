# [[file:../../org/pyretechnics.org::py-types-py][py-types-py]]
import cython as cy

#==============================================================
# Runtime-defined type aliases
#==============================================================

pyidx     = cy.typedef(cy.Py_ssize_t)
vec_xy    = cy.typedef(tuple[cy.float, cy.float])
vec_xyz   = cy.typedef(tuple[cy.float, cy.float, cy.float])
coord_yx  = cy.typedef(tuple[pyidx, pyidx])
coord_tyx = cy.typedef(tuple[pyidx, pyidx, pyidx])
# py-types-py ends here
