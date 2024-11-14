# [[file:../../org/pyretechnics.org::types-py][types-py]]
import cython
import cython as cy

#==============================================================
# Runtime-defined type aliases
#==============================================================

if cy.compiled:
    from cython.cimports.pyretechnics.types import pyidx, vec_xy, vec_xyz
else:
    pyidx   = cy.typedef(cy.Py_ssize_t)
    vec_xy  = cy.typedef(tuple[cy.float, cy.float])
    vec_xyz = cy.typedef(tuple[cy.float, cy.float, cy.float])
# types-py ends here
