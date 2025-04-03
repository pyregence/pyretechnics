import cython
import cython as cy
import numpy as np
from numpy.random import Generator
if cython.compiled:
    from cython.cimports.pyretechnics.cy_types import pyidx
else:
    from pyretechnics.py_types import pyidx


@cy.cclass
class BufferedRandGen:
    numpy_rand     : Generator
    uniform_buf    : cy.float[::1]
    uniform_pos    : pyidx
    normal_buf     : cy.float[::1]
    normal_pos     : pyidx
    poisson16_buf  : cy.longlong[::1] # An array of Poisson(16) draws.
    poisson16_pos  : pyidx
    poisson1_buf   : cy.longlong[::1] # An array of Poisson(1) draws.
    poisson1_pos   : pyidx
    poisson_exp_buf: cy.double[::1]   # An array of Exponential(1) draws, with the first value potentially modified.
    poisson_exp_pos: pyidx
    # NOTE: INVARIANT poisson_exp_buf[poisson_exp_pos] is the distance to the next Poisson Process point.


    def __cinit__(self, numpy_rand):
        self.numpy_rand = numpy_rand
        __reset_uniform_buffer(self)
        __reset_normal_buffer(self)
        self.poisson16_buf = self.numpy_rand.poisson(lam=16.0, size=1024).astype(np.int64)
        self.poisson1_buf  = self.numpy_rand.poisson(lam=1.0, size=1024).astype(np.int64)
        __reset_exp_buffer(self)


    @cy.cfunc
    @cy.exceptval(check=False)
    def next_poisson(self, M: cy.double) -> cy.longlong:
        """
        A method for efficiently drawing from a Poisson distribution of (very) small mean,
        doing fewer random draws than method calls.

        Returns a random draw from a Poisson distribution of mean M.
        """
        if M > 0.0:
            ret: cy.longlong = 0
            # Theorem: A sum of independent Poisson variables is Poisson-distributed,
            #          and its mean is the sum of the means.
            # Draw repeatedly from Poisson(16)
            while M >= 16.0:
                ret += __next_poisson16(self)
                M   -= 16.0
            # Draw repeatedly from Poisson(1)
            while M >= 1.0:
                ret += __next_poisson1(self)
                M   -= 1.0
            # Now we draw efficiently from a Poisson distribution of fractional mean.
            # This relies on using a sequence of arrival times distributed i.i.d. as Exponential(1).
            ret_frac: cy.longlong = 0
            pos     : pyidx       = self.poisson_exp_pos
            next_Sk : cy.double   = self.poisson_exp_buf[pos]
            while M >= next_Sk: # This will be rare.
                ret_frac += 1
                M        -= next_Sk
                pos      += 1
                if pos >= 64:
                    __reset_exp_buffer(self)
                    pos = self.poisson_exp_pos
                next_Sk = self.poisson_exp_buf[pos]
            self.poisson_exp_buf[pos] = next_Sk - M
            if ret_frac > 0: # We moved the cursor.
                self.poisson_exp_pos = pos
            ret += ret_frac
            return ret
        else:
            return 0


    @cy.cfunc
    @cy.exceptval(check=False)
    def next_uniform(self) -> cy.float:
        if not(self.uniform_pos < 1024):
            __reset_uniform_buffer(self)
        ret: cy.float     = self.uniform_buf[self.uniform_pos]
        self.uniform_pos += 1
        return ret


    @cy.cfunc
    @cy.exceptval(check=False)
    def next_normal(self) -> cy.float:
        if not(self.normal_pos < 1024):
            __reset_normal_buffer(self)
        ret: cy.float    = self.normal_buf[self.normal_pos]
        self.normal_pos += 1
        return ret


@cy.cfunc
@cy.exceptval(check=False)
def __reset_exp_buffer(self: BufferedRandGen) -> cy.void:
    self.poisson_exp_buf = self.numpy_rand.exponential(size=64)
    self.poisson_exp_pos = 0


@cy.cfunc
@cy.exceptval(check=False)
def __reset_uniform_buffer(self: BufferedRandGen) -> cy.void:
    self.uniform_buf = self.numpy_rand.uniform(size=1024).astype(np.float32)
    self.uniform_pos = 0


@cy.cfunc
@cy.exceptval(check=False)
def __reset_normal_buffer(self: BufferedRandGen) -> cy.void:
    self.normal_buf = self.numpy_rand.normal(size=1024).astype(np.float32)
    self.normal_pos = 0


@cy.cfunc
@cy.exceptval(check=False)
def __next_poisson1(self: BufferedRandGen) -> cy.longlong:
    if not(self.poisson1_pos < 1024):
        self.poisson1_buf = self.numpy_rand.poisson(lam=1.0, size=1024).astype(np.int64)
        self.poisson1_pos = 0
    ret: cy.longlong   = self.poisson1_buf[self.poisson1_pos]
    self.poisson1_pos += 1
    return ret


@cy.cfunc
@cy.exceptval(check=False)
def __next_poisson16(self: BufferedRandGen) -> cy.longlong:
    if not(self.poisson16_pos < 1024):
        self.poisson16_buf = self.numpy_rand.poisson(lam=16.0, size=1024).astype(np.int64)
        self.poisson16_pos = 0
    ret: cy.longlong    = self.poisson16_buf[self.poisson16_pos]
    self.poisson16_pos += 1
    return ret
