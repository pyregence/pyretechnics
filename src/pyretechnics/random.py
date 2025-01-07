import cython
if cython.compiled:
    from cython.cimports.pyretechnics.math import log
    from cython.cimports.pyretechnics.cy_types import pyidx
else:
    from pyretechnics.py_types import pyidx
    from cython.cimports.pyretechnics.math import log

import cython as cy
import numpy as np
from numpy.random import Generator

@cy.cclass
class BufferedRandGen:
    numpy_rand: Generator
    
    uniform_pos: pyidx
    uniform_buf: cy.float[:]

    normal_pos: pyidx
    normal_buf: cy.float[:]

    def __cinit__(self, numpy_rand):
        self.numpy_rand = numpy_rand
        _reset_uniform_buffer(self)
        _reset_normal_buffer(self)

    @cy.ccall
    def next_uniform(self) -> cy.float:
        if not(self.uniform_pos < 1024):
            _reset_uniform_buffer(self)
        ret: cy.float = self.uniform_buf[self.uniform_pos]
        self.uniform_pos += 1
        return ret

    @cy.ccall
    def next_normal(self) -> cy.float:
        if not(self.normal_pos < 1024):
            _reset_normal_buffer(self)
        ret: cy.float = self.normal_buf[self.normal_pos]
        self.normal_pos += 1
        return ret     


            
@cy.ccall
def _reset_uniform_buffer(self: BufferedRandGen) -> cy.void:
    bfr: cy.float[:] = self.numpy_rand.uniform(size=1024).astype(np.float32)
    self.uniform_buf = bfr
    self.uniform_pos = 0


@cy.ccall
def _reset_normal_buffer(self: BufferedRandGen) -> cy.void:
    bfr: cy.float[:] = self.numpy_rand.normal(size=1024).astype(np.float32)
    self.normal_buf = bfr
    self.normal_pos = 0

