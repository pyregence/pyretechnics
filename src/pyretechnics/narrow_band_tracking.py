"""
Custom data structures for keeping track of cells in the narrow band.

The narrow band is defined as the set of _tracked_ cells:
- a cell is _tracked_ iff it lies within a square buffer around a _frontier cell_;
- a _frontier cell_ is a cell that has a neighbor of opposite phi sign,
i.e. both cells lie on opposite sides of the fire front.

`NarrowBandTracker` maintains the set of tracked cells by mapping each track cell to a count,
the count of the frontier cells related to it.
The intended use is that the surrounding algorithm will follow the additions and deletions
to the set of frontier cells, and correspondingly call the functions
`incr_square_around()` and `decr_square_around()`.
`TrackedCellsIterator` allows to iterate over tracked cells
(more efficiently than a regular Python iterator).
"""

import cython
if cython.compiled:
    from cython.cimports.pyretechnics.cy_types import pyidx, coord_yx
else:
    from pyretechnics.py_types import pyidx, coord_yx
import cython as cy
import sortedcontainers as sortc    


@cy.cclass
class CellsCountSegment:
    """
    For internal use of NarrowBandTracker.

    A CellsCountSegment maps each (y, x) for x from x0 (included) to x0+16 (excluded)
    to a positive integer.
    """

    y: cy.int
    x0: cy.int
    counts: cy.ushort[16]


@cy.cfunc
@cy.exceptval(check=False)
@cy.boundscheck(False)
@cy.wraparound(False)
@cy.inline
def segm_is_empty(segm: CellsCountSegment) -> cy.bint:
    for k in range(16):
        if segm.counts[k] > 0:
            return False
    return True


@cy.cfunc
@cy.exceptval(check=False)
@cy.boundscheck(False)
@cy.wraparound(False)
@cy.inline
def segm_is_pos_at(segm: CellsCountSegment, k: pyidx) -> cy.bint:
    return (segm.counts[k] > 0)


@cy.ccall
def make_CellsCountSegment(y: cy.int, x0: cy.int) -> CellsCountSegment:
    """
    Creates an empty CellsCountSegment (with all counts to zero).
    """
    ret: CellsCountSegment = CellsCountSegment()
    ret.y = y
    ret.x0 = x0
    return ret


@cy.ccall
@cy.exceptval(-1)
def _find_first_k(segm: CellsCountSegment) -> cy.int:
    k: cy.int
    for k in range(16):
        if segm_is_pos_at(segm, k):
            return k
    raise ValueError("Segment is empty!")



@cy.cclass
class NarrowBandTracker:
    """
    A custom data structure for keeping track of cells in the narrow band efficiently.
    Logically, this class acts a sorted dict from (y, x) -> n,
    in which n is the number of frontier cells which buffer contains (y, x).
    """
    n_tracked_cells: pyidx

    y_high: pyidx
    x_high: pyidx
    # INVARIANT for all (y, x) in the keys, ys_offset <= y < ys_offset + len(ys_list)
    ys_offset: cy.int
    # OPTIM maybe we can make this faster by making ys_list a raw array.
    ys_list: list # ys_list[i] logically represents the row for y := ys_offset + i.
    # INVARIANT ys_list[i] is either None or a non-empty SortedDict s;
    # ys_list[i] is None iff the row is empty.
    # INVARIANT the SortecDict s := ys_list[i] maps each (x // 16) value to a non-empty CellsCountSegment segm,
    # i.e. one of the counts in segm.counts will be > 0.

    _rows_count: cy.int # INVARIANT the count of y_idx for which ys_list[y_idx] is not None. Therefore the count of distinct y values with positive counts.


@cy.ccall
def new_NarrowBandTracker(y_high: pyidx, x_high: pyidx) -> NarrowBandTracker:
    ret: NarrowBandTracker = NarrowBandTracker()
    ret.n_tracked_cells = 0
    ret.y_high = y_high
    ret.x_high = x_high
    return ret


# NOTE writing these are like methods, but writing them as functions keeps them private and simplifies pxd declarations.
@cy.ccall
def _ensure_y(tracked_cells: NarrowBandTracker, y: cy.int):
    """
    Ensures that index y is covered by ys_offset and ys_list.
    """
    if tracked_cells.ys_list is None:
        tracked_cells.ys_offset = y
        tracked_cells.ys_list = [None] * 8
    n: cy.int = len(tracked_cells.ys_list)
    while y < tracked_cells.ys_offset: # Double the size, extending left
        ys_list_new: list = [None] * (2*n)
        ys_list_new[n:(2*n)] = tracked_cells.ys_list
        tracked_cells.ys_offset -= n
        n = 2*n
        tracked_cells.ys_list = ys_list_new
    while y >= tracked_cells.ys_offset + n: # Double the size, extending right
        tracked_cells.ys_list[n:(2*n)] = [None] * n
        n = 2*n



# OPTIM list and dict lookup operations require us to allocate tons of boxed integers on the heap.
# We can maybe make this faster by caching those integers in a native array.
# That said, I expect Python to keep doing reference-counting mutations during function calls...

@cy.ccall
def incr_y_segment(tracked_cells: NarrowBandTracker, y: cy.int, x_start: cy.int, segm_length: cy.int):
    _ensure_y(tracked_cells, y)
    y_idx: object = y - tracked_cells.ys_offset
    s: sortc.SortedDict = tracked_cells.ys_list[y_idx]
    if s is None:
        s = sortc.SortedDict()
        tracked_cells.ys_list[y_idx] = s
        tracked_cells._rows_count += 1
    xk: cy.int = (x_start // 16) - 1
    segm: object = None
    n_added: pyidx = 0
    i: cy.int
    for i in range(segm_length):
        if (x_start + i)//16 != xk:
            xk = (x_start + i)//16
            segm = s.get(xk)
            if segm is None:
                segm = make_CellsCountSegment(y, xk*16)
                s[xk] = segm
        segment: CellsCountSegment = segm
        k: pyidx = (x_start + i) % 16
        old_count: cy.ushort = segment.counts[k]
        segment.counts[k] = old_count + 1
        n_added += (old_count == 0)
    tracked_cells.n_tracked_cells += n_added

@cy.ccall
def decr_y_segment(tracked_cells: NarrowBandTracker, y: cy.int, x_start: cy.int, segm_length: cy.int):
    y_idx: object = y - tracked_cells.ys_offset
    s: object = tracked_cells.ys_list[y_idx]
    if s is None:
        raise ValueError("Empty SortedDict")
    xk: cy.int = (x_start // 16)
    segm: cy.Optional[CellsCountSegment] = s.get(xk)
    if segm is None:
        raise ValueError("No Segment found!")
    n_removed: pyidx = 0
    i: cy.int
    for i in range(segm_length):
        k: cy.int = (x_start + i) % 16
        segment: CellsCountSegment = segm
        old_count: cy.ushort = segment.counts[k]
        segment.counts[k] = old_count - 1
        n_removed += (old_count == 1)
        if segm_is_empty(segment):
            del s[xk]
        if (x_start + i + 1)//16 != xk: # we're about to change segment
            xk = (x_start + i + 1)//16
            segm = s.get(xk)
    if len(s) == 0:
        tracked_cells.ys_list[y_idx] = None
        tracked_cells._rows_count -= 1
    tracked_cells.n_tracked_cells -= n_removed

@cy.ccall
@cy.exceptval(check=False)
def resolve_truncated_x_segment(tracked_cells: NarrowBandTracker, x: cy.int, buffer_width: cy.int) -> cy.tuple[cy.int, cy.int]:
    x_start: cy.int = max(0, x - buffer_width)
    x_end: cy.int = min(tracked_cells.x_high - 1, x + buffer_width)
    width_truncated: cy.int = x_end - x_start
    return (x_start, width_truncated)


@cy.ccall
def incr_square_around(tracked_cells: NarrowBandTracker, y: cy.int, x: cy.int, buffer_width: cy.int):
    """
    Mutates tracked_cells by incrementing all cells within `buffer_width` of (y, x).

    For example, if buffer_width = 3, then 7*7=49 counts are updated.
    """
    width: cy.int = 2*buffer_width + 1
    i: cy.int
    x_start, width_truncated = resolve_truncated_x_segment(tracked_cells, x, buffer_width)
    for i in range(width):
        y1: cy.int = y - buffer_width + i
        if (y1 >= 0) and (y1 < tracked_cells.y_high):
            incr_y_segment(tracked_cells, y1, x_start, width_truncated)


@cy.ccall
def decr_square_around(tracked_cells: NarrowBandTracker, y: cy.int, x: cy.int, buffer_width: cy.int):
    """
    Mutates tracked_cells by decrementing all cells within `buffer_width` of (y, x).

    If buffer_width = 3, then 7*7=49 counts are updated.

    The mutated cells MUST all have a positive count before this function is called.
    """
    width: cy.int = 2*buffer_width + 1
    i: cy.int
    x_start, width_truncated = resolve_truncated_x_segment(tracked_cells, x, buffer_width)
    for i in range(width):
        y1: cy.int = y - buffer_width + i
        if (y1 >= 0) and (y1 < tracked_cells.y_high):
            decr_y_segment(tracked_cells, y1, x_start, width_truncated)


@cy.ccall
def nonempty_tracked_cells(tracked_cells: NarrowBandTracker) -> cy.bint:
    return tracked_cells.n_tracked_cells > 0


# Iterating over cells
# NOTE @cy.cclass is not applicable here, due to yield.
def iterate_segments(tracked_cells: NarrowBandTracker):
    """
    Returns a Python iterator over the CellsCountSegment-s in `tracked_cells`.
    """
    ys_list: list = tracked_cells.ys_list
    if ys_list is not None:
        for s in ys_list:
            if s is not None:
                for segm in s.values():
                    yield segm


@cy.cclass
class TrackedCellsIterator:
    """
    An "Iterator-pattern" object for iterating efficiently over the cells of a NarrowBandTracker.
    """
    segm_iter: object # A Python Iterator of CellsCountSegment. INVARIANT None iff iteration is exhausted.
    current_segm: CellsCountSegment # INVARIANT never None
    current_k: cy.int # INVARIANT points to the index for the next return cell.

    # NOTE this implementation relies on the fact that a CellsCountSegment is never empty.
    def __cinit__(self, segm_iter):
        self.segm_iter = segm_iter
        try:
            self.current_segm = next(segm_iter) # OPTIM maybe we could make even this faster by iterating directly in NarrowBandTracker instead of relying on a Python iterator, effectively reimplementing iterate_segments() without yield. Basic testing suggests of 4x potential efficiency gain by getting rid of yield.
            self.current_k = _find_first_k(self.current_segm)
        except StopIteration: # Only if the iterator started empty!
            self.segm_iter = None
            self.current_segm = make_CellsCountSegment(-1, -1)
            self.current_k = 0
        

    @cy.ccall
    def has_next(self) -> cy.bint: # OPTIM it may be faster to maintain an iteration index and compare it to n_tracked_cells.
        return (self.segm_iter is not None)

    @cy.ccall
    def next_cell(self) -> coord_yx:
        segm: CellsCountSegment = self.current_segm
        y: pyidx = segm.y
        x: pyidx = segm.x0 + self.current_k
        ret: coord_yx = (y, x)
        _resolve_next(self)
        return ret


@cy.ccall
def _resolve_next(tcitr: TrackedCellsIterator) -> cy.bint:
    # At the end of this, either segm_iter is None (exhausted),
    # or current_segm and current_k point to the next cell.
    # Returns whether iteration exhausted.
    current_segm: CellsCountSegment = tcitr.current_segm
    k: cy.int = tcitr.current_k + 1
    while k < 16:
        if segm_is_pos_at(current_segm, k):
            tcitr.current_k = k
            return False
        k += 1
    # We exhausted the segment, beginning the next one:
    try:
        tcitr.current_segm = next(tcitr.segm_iter)
    except StopIteration:
        tcitr.segm_iter = None
        return True
    tcitr.current_k = _find_first_k(tcitr.current_segm)
    return False


@cy.ccall
def tracked_cells_iterator(tracked_cells: NarrowBandTracker) -> TrackedCellsIterator:
    return TrackedCellsIterator(iterate_segments(tracked_cells))
