from pyretechnics.cy_types cimport pyidx, coord_yx

cdef class CellsCountSegment:
    cdef int y
    cdef int x0
    cdef unsigned short[16] counts

    cpdef bint is_pos_at(self, pyidx k)
    cpdef bint is_empty(self)


cdef class NarrowBandTracker:
    cdef pyidx y_high
    cdef pyidx x_high
    cdef int ys_offset
    cdef list ys_list
    cdef int _rows_count

cpdef NarrowBandTracker new_NarrowBandTracker(pyidx y_high, pyidx x_high)


cdef class TrackedCellsIterator:
    cdef object segm_iter
    cdef CellsCountSegment current_segm
    cdef int current_k
    # def __cinit__(self, segm_iter)
    cpdef bint has_next(self)
    cpdef coord_yx next_cell(self)


cpdef TrackedCellsIterator tracked_cells_iterator(NarrowBandTracker tracked_cells)
cpdef bint nonempty_tracked_cells(NarrowBandTracker tracked_cells)


cpdef incr_square_around(NarrowBandTracker tracked_cells, int y, int x, int buffer_width)
cpdef decr_square_around(NarrowBandTracker tracked_cells, int y, int x, int buffer_width)

