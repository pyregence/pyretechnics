# [[file:../../org/pyretechnics.org::types-py][types-py]]
import cython
import cython as cy

#==============================================================
# Runtime-defined type aliases
#==============================================================

if cy.compiled:
    from cython.cimports.pyretechnics.types import pyidx, vec_xy, vec_xyz, coord_yx, coord_tyx
else:
    # TODO: Maybe pyidx should be cy.int?
    # TODO: Maybe we should use C arrays instead of tuples(structs)?
    pyidx     = cy.typedef(cy.Py_ssize_t)
    vec_xy    = cy.typedef(tuple[cy.float, cy.float])
    vec_xyz   = cy.typedef(tuple[cy.float, cy.float, cy.float])
    coord_yx  = cy.typedef(tuple[pyidx, pyidx])
    coord_tyx = cy.typedef(tuple[pyidx, pyidx, pyidx])
# types-py ends here
