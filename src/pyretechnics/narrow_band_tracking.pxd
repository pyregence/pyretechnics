from pyretechnics.cy_types cimport pyidx, coord_yx

cdef class CellsCountSegment:
    cdef int y
    cdef int x0
    cdef unsigned short[16] counts

cdef bint segment_is_empty(CellsCountSegment segment) noexcept
cdef bint segment_is_pos_at(CellsCountSegment segment, pyidx k) noexcept
cdef CellsCountSegment make_CellsCountSegment(pyidx y, pyidx x0)

cdef class NarrowBandTracker:
    cdef public pyidx n_tracked_cells # TODO: Why is this marked public?
    cdef pyidx y_high
    cdef pyidx x_high
    cdef pyidx ys_offset
    cdef list ys_list
    cdef pyidx _rows_count

cdef NarrowBandTracker new_NarrowBandTracker(pyidx y_high, pyidx x_high)
cdef void inc_y_segment(NarrowBandTracker tracked_cells, pyidx y, pyidx x_start, int segment_length) noexcept
cdef void dec_y_segment(NarrowBandTracker tracked_cells, pyidx y, pyidx x_start, int segment_length) noexcept
cdef (pyidx, pyidx) resolve_truncated_x_segment(NarrowBandTracker tracked_cells, pyidx x, pyidx buffer_width) noexcept
cdef void inc_square_around(NarrowBandTracker tracked_cells, pyidx y, pyidx x, pyidx buffer_width) noexcept
cdef void dec_square_around(NarrowBandTracker tracked_cells, pyidx y, pyidx x, pyidx buffer_width) noexcept
cdef bint nonempty_tracked_cells(NarrowBandTracker tracked_cells) noexcept

cdef class TrackedCellsIterator:
    cdef object segment_iter
    cdef CellsCountSegment current_segment
    cdef pyidx current_k
    cdef bint has_next(TrackedCellsIterator self) noexcept
    cdef coord_yx next_cell(TrackedCellsIterator self) noexcept

cdef TrackedCellsIterator tracked_cells_iterator(NarrowBandTracker tracked_cells)
