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
`inc_square_around()` and `dec_square_around()`.
`TrackedCellsIterator` allows to iterate over tracked cells
(more efficiently than a regular Python iterator).
"""
import cython
import cython as cy
from sortedcontainers import SortedDict
if cython.compiled:
    from cython.cimports.pyretechnics.cy_types import pyidx, coord_yx
else:
    from pyretechnics.py_types import pyidx, coord_yx

#===================================================================
# Class CellsCountSegment
#===================================================================

@cy.cclass
class CellsCountSegment:
    """
    For internal use by the NarrowBandTracker.

    A CellsCountSegment maps each (y, x) for x from x0 (included) to x0+16 (excluded)
    to a positive integer.
    """
    y     : pyidx
    x0    : pyidx
    counts: cy.ushort[16]


@cy.cfunc
@cy.exceptval(check=False)
def segment_is_empty(segment: CellsCountSegment) -> cy.bint:
    counts: cy.ushort[16] = segment.counts
    k: pyidx
    for k in range(16):
        if counts[k] > 0:
            return False
    return True


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def segment_is_pos_at(segment: CellsCountSegment, k: pyidx) -> cy.bint:
    return segment.counts[k] > 0


# FIXME: Why isn't this a class __init__ constructor? (Val: ensuring private functions, simplifying pxd declarations)
@cy.cfunc
def make_CellsCountSegment(y: pyidx, x0: pyidx) -> CellsCountSegment:
    """
    Creates an empty CellsCountSegment (with all counts set to zero).
    """
    ret: CellsCountSegment = CellsCountSegment()
    ret.y  = y
    ret.x0 = x0
    return ret


@cy.cfunc
def __find_first_k(segment: CellsCountSegment) -> pyidx:
    k: pyidx
    for k in range(16):
        if segment_is_pos_at(segment, k):
            return k
    raise ValueError("Segment is empty!")

#===================================================================
# Class NarrowBandTracker
#===================================================================

# TODO: OPTIM Maybe we can make this faster by making ys_list a raw array.
@cy.cclass
class NarrowBandTracker:
    """
    A custom data structure for keeping track of cells in the narrow band efficiently.
    Logically, this class acts as a sorted dict from (y, x) -> n,
    in which n is the number of frontier cells whose buffer contains (y, x).

    Invariants:
    - for all (y, x) in the keys, ys_offset <= y < ys_offset + len(ys_list)
    - ys_list[i] logically represents the row for y := ys_offset + i
    - ys_list[i] is either None or a non-empty SortedDict s
    - ys_list[i] is None iff the row is empty
    - the SortedDict s := ys_list[i] maps each (x // 16) value to a non-empty CellsCountSegment segment,
      i.e. one of the counts in segment.counts will be > 0
    - _rows_count is the number of entries in ys_list containing non-empty SortedDicts
    """
    num_tracked_cells: cy.int
    y_high           : pyidx
    x_high           : pyidx
    ys_offset        : pyidx
    ys_list          : list
    _rows_count      : cy.int


# FIXME: Why isn't this a class __init__ constructor? (Val: ensuring private functions, simplifying pxd declarations)
@cy.cfunc
def new_NarrowBandTracker(y_high: pyidx, x_high: pyidx) -> NarrowBandTracker:
    ret: NarrowBandTracker = NarrowBandTracker()
    ret.num_tracked_cells  = 0
    ret.y_high             = y_high
    ret.x_high             = x_high
    return ret


@cy.cfunc
@cy.exceptval(check=False)
def __ensure_y(tracked_cells: NarrowBandTracker, y: pyidx) -> cy.void:
    """
    Ensures that index y is covered by ys_offset and ys_list.
    """
    # Ensure that the ys_offset and ys_list fields of tracked_cells are initialized
    if tracked_cells.ys_list is None:
        tracked_cells.ys_offset = y
        tracked_cells.ys_list   = [None] * 8
    # Grow the tracked_cells contents to the left or right if necessary
    n: cy.int = len(tracked_cells.ys_list)
    while y < tracked_cells.ys_offset:
        # Double the size, extending left
        tracked_cells.ys_list    = ([None] * n) + tracked_cells.ys_list
        tracked_cells.ys_offset -= n # FIXME: What happens when this becomes negative?
        n = 2*n
    while y >= tracked_cells.ys_offset + n:
        # Double the size, extending right
        tracked_cells.ys_list = tracked_cells.ys_list + ([None] * n)
        n = 2*n


# TODO: OPTIM list and dict lookup operations require us to allocate
#       tons of boxed integers on the heap. We can maybe make this
#       faster by caching those integers in a native array. That said,
#       I expect Python to keep doing reference-counting mutations
#       during function calls...
@cy.cfunc
@cy.exceptval(check=False)
def inc_y_segment(tracked_cells: NarrowBandTracker, y: pyidx, x_start: pyidx, segment_length: cy.int) -> cy.void:
    __ensure_y(tracked_cells, y)
    y_idx: int        = y - tracked_cells.ys_offset
    s    : SortedDict = tracked_cells.ys_list[y_idx]
    if s is None:
        s                             = SortedDict()
        tracked_cells.ys_list[y_idx]  = s
        tracked_cells._rows_count    += 1
    xk       : pyidx                  = (x_start // 16) - 1
    segment0 : CellsCountSegment|None = None
    segment  : CellsCountSegment
    n_added  : cy.int                 = 0
    old_count: cy.ushort
    k        : pyidx
    i        : pyidx
    for i in range(segment_length):
        if (x_start + i) // 16 != xk:
            xk       = (x_start + i) // 16
            segment0 = s.get(xk)
            if segment0 is None:
                segment0 = make_CellsCountSegment(y, xk*16)
                s[xk]    = segment0
        segment            = segment0
        k                  = (x_start + i) % 16
        old_count          = segment.counts[k]
        segment.counts[k]  = old_count + 1
        n_added           += (old_count == 0)
    tracked_cells.num_tracked_cells += n_added


@cy.cfunc
def dec_y_segment(tracked_cells: NarrowBandTracker, y: pyidx, x_start: pyidx, segment_length: cy.int) -> cy.void:
    y_idx: int        = y - tracked_cells.ys_offset
    s    : SortedDict = tracked_cells.ys_list[y_idx]
    if s is None:
        raise ValueError("Empty SortedDict!")
    xk      : pyidx                  = (x_start // 16)
    segment0: CellsCountSegment|None = s.get(xk)
    if segment0 is None:
        raise ValueError("No CellsCountSegment found!")
    segment  : CellsCountSegment
    n_removed: cy.int = 0
    old_count: cy.ushort
    k        : pyidx
    i        : pyidx
    for i in range(segment_length):
        k                  = (x_start + i) % 16
        segment            = segment0
        old_count          = segment.counts[k]
        segment.counts[k]  = old_count - 1
        n_removed         += (old_count == 1)
        if segment_is_empty(segment):
            del s[xk]
        if (x_start + i + 1) // 16 != xk:
            xk       = (x_start + i + 1) // 16
            segment0 = s.get(xk)
    if len(s) == 0:
        tracked_cells.ys_list[y_idx]  = None
        tracked_cells._rows_count    -= 1
    tracked_cells.num_tracked_cells  -= n_removed


@cy.cfunc
@cy.exceptval(check=False)
def resolve_truncated_x_segment(tracked_cells: NarrowBandTracker,
                                x            : pyidx,
                                buffer_width : pyidx) -> tuple[pyidx, pyidx]:
    x_start        : pyidx = max(0, x - buffer_width)
    x_end          : pyidx = min(tracked_cells.x_high - 1, x + buffer_width)
    width_truncated: pyidx = x_end - x_start
    return (x_start, width_truncated)


@cy.cfunc
@cy.exceptval(check=False)
def inc_square_around(tracked_cells: NarrowBandTracker, y: pyidx, x: pyidx, buffer_width: pyidx) -> cy.void:
    """
    Mutates tracked_cells by incrementing all cells within `buffer_width` of (y, x).

    For example, if buffer_width = 3, then 7*7=49 counts are updated.
    """
    (x_start, width_truncated) = resolve_truncated_x_segment(tracked_cells, x, buffer_width)
    width: pyidx               = 2 * buffer_width + 1
    i    : pyidx
    for i in range(width):
        y1: pyidx = y - buffer_width + i
        if (y1 >= 0) and (y1 < tracked_cells.y_high):
            inc_y_segment(tracked_cells, y1, x_start, width_truncated)


@cy.cfunc
@cy.exceptval(check=False)
def dec_square_around(tracked_cells: NarrowBandTracker, y: pyidx, x: pyidx, buffer_width: pyidx) -> cy.void:
    """
    Mutates tracked_cells by decrementing all cells within `buffer_width` of (y, x).

    If buffer_width = 3, then 7*7=49 counts are updated.

    The mutated cells MUST all have a positive count before this function is called.
    """
    (x_start, width_truncated) = resolve_truncated_x_segment(tracked_cells, x, buffer_width)
    width: pyidx               = 2 * buffer_width + 1
    i    : pyidx
    for i in range(width):
        y1: pyidx = y - buffer_width + i
        if (y1 >= 0) and (y1 < tracked_cells.y_high):
            dec_y_segment(tracked_cells, y1, x_start, width_truncated)


@cy.cfunc
@cy.inline
@cy.exceptval(check=False)
def nonempty_tracked_cells(tracked_cells: NarrowBandTracker) -> cy.bint:
    return tracked_cells.num_tracked_cells > 0


# Iterating over cells
# NOTE: @cy.cfunc is not applicable here, due to yield.
def iterate_segments(tracked_cells: NarrowBandTracker):
    """
    Returns a Python iterator over the CellsCountSegment-s in `tracked_cells`.
    """
    ys_list: list = tracked_cells.ys_list
    if ys_list is not None:
        for s in ys_list:
            if s is not None:
                for segment in s.values():
                    yield segment


# TODO: OPTIM Maybe we could make this even faster by iterating
#       directly in NarrowBandTracker instead of relying on a Python
#       iterator, effectively reimplementing iterate_segments()
#       without yield. Basic testing suggests a 4x potential
#       efficiency gain by getting rid of yield!
@cy.cclass
class TrackedCellsIterator:
    """
    An "Iterator-pattern" object for iterating efficiently over the cells of a NarrowBandTracker.
    """
    segment_iter   : object # A Python Iterator of CellsCountSegment. INVARIANT None iff iteration is exhausted.
    current_segment: CellsCountSegment # INVARIANT never None
    current_k      : pyidx # INVARIANT points to the index for the next cell to return.


    # NOTE: This implementation relies on the fact that a CellsCountSegment is never empty.
    def __cinit__(self, segment_iter):
        self.segment_iter = segment_iter
        try:
            self.current_segment = next(segment_iter)
            self.current_k       = __find_first_k(self.current_segment)
        except StopIteration: # Only if the iterator started empty!
            self.segment_iter    = None
            self.current_segment = make_CellsCountSegment(-1, -1)
            self.current_k       = 0


    # TODO: OPTIM It may be faster to maintain an iteration index and compare it to num_tracked_cells.
    @cy.cfunc
    @cy.exceptval(check=False)
    def has_next(self) -> cy.bint:
        return self.segment_iter is not None


    # FIXME: Cython uninitialized warning
    @cy.cfunc
    @cy.exceptval(check=False)
    def next_cell(self) -> coord_yx:
        segment: CellsCountSegment = self.current_segment
        y      : pyidx             = segment.y
        x      : pyidx             = segment.x0 + self.current_k
        ret    : coord_yx          = (y, x)
        __resolve_next(self)
        return ret


@cy.cfunc
@cy.exceptval(check=False)
def __resolve_next(tc_iter: TrackedCellsIterator) -> cy.bint:
    # At the end of this, either segment_iter is None (exhausted),
    # or current_segment and current_k point to the next cell.
    # Returns True if iteration was exhausted and False otherwise.
    current_segment: CellsCountSegment = tc_iter.current_segment
    k              : pyidx             = tc_iter.current_k + 1
    while k < 16:
        if segment_is_pos_at(current_segment, k):
            tc_iter.current_k = k
            return False
        k += 1
    # We exhausted the segment, beginning the next one:
    try:
        tc_iter.current_segment = next(tc_iter.segment_iter)
    except StopIteration:
        tc_iter.segment_iter = None
        return True
    tc_iter.current_k = __find_first_k(tc_iter.current_segment)
    return False


@cy.cfunc
@cy.inline
def tracked_cells_iterator(tracked_cells: NarrowBandTracker) -> TrackedCellsIterator:
    return TrackedCellsIterator(iterate_segments(tracked_cells))
