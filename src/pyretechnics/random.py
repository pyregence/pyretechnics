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

    poisson16buf: cy.long[:] # An array of Poisson(16) draws.
    poisson16pos: pyidx
    poisson1buf: cy.long[:] # An array of Poisson(1) draws.
    poisson1pos: pyidx
    poisson_exp_buf: cy.double[:] # An array of Exponential(1) draws, with the first value potentially modified.
    poisson_exp_pos: pyidx # Cursor in the above array
    # INVARIANT poisson_exp_buf[poisson_exp_pos] is the distance to the next Poisson Process point.

    def __cinit__(self, numpy_rand):
        self.numpy_rand = numpy_rand
        _reset_uniform_buffer(self)
        _reset_normal_buffer(self)
        self.poisson16buf = self.numpy_rand.poisson(lam=16., size=1024)
        self.poisson1buf = self.numpy_rand.poisson(lam=1., size=1024)
        _reset_exp_buffer(self)
        

    @cy.ccall
    def next_poisson(self, M: cy.double) -> pyidx:
        """
        A method for efficiently drawing from a Poisson distribution of (very) small mean,
        doing fewer random draws than method calls.

        Returns a random draw from a Poisson distribution of mean M.
        """
        if M > 0:
            ret: pyidx = 0
            # Theorem: a sum of independent Poisson variables is Poisson-distributed, and its mean is the sum of the means.
            # Draw repeatedly from Poisson(16)
            while M >= 16:
                ret += _next_poisson16(self)
                M -= 16
            # Draw repeatedly from Poisson(1)
            while M >= 1:
                ret += _next_poisson1(self)
                M -= 1
            # Now we draw efficiently from a Poisson distribution of fractional mean.
            # This relies on using a sequence of arrival times distributed i.i.d. as Exponential(1).
            ret_frac: pyidx = 0
            pos: pyidx = self.poisson_exp_pos
            next_Sk: cy.double = self.poisson_exp_buf[pos]
            while M >= next_Sk: # This will be rare.
                ret_frac += 1
                M -= next_Sk
                pos += 1
                if pos >= 64:
                    _reset_exp_buffer(self)
                    pos = self.poisson_exp_pos
                next_Sk = self.poisson_exp_buf[pos]
            self.poisson_exp_buf[pos] = next_Sk - M
            if ret_frac > 0: # we moved the cursor
                self.poisson_exp_pos = pos
            ret += ret_frac
            return ret
        else:
            return 0


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
def _reset_exp_buffer(self: BufferedRandGen) -> cy.void:
    bfr: cy.double[:] = self.numpy_rand.exponential(size=64)
    self.poisson_exp_buf = bfr
    self.poisson_exp_pos = 0


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


@cy.ccall
def _next_poisson1(self: BufferedRandGen) -> pyidx:
    if not(self.poisson1pos < 1024):
        bfr: cy.long[:] = self.numpy_rand.poisson(lam=1., size=1024)
        self.poisson1buf = bfr
        self.poisson1pos = 0
    ret: pyidx = self.poisson1buf[self.poisson1pos]
    self.poisson1pos += 1
    return ret

@cy.ccall
def _next_poisson16(self: BufferedRandGen) -> pyidx:
    if not(self.poisson16pos < 1024):
        bfr: cy.long[:] = self.numpy_rand.poisson(lam=16., size=1024)
        self.poisson16buf = bfr
        self.poisson16pos = 0
    ret: pyidx = self.poisson16buf[self.poisson16pos]
    self.poisson16pos += 1
    return ret
