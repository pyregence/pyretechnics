# [[file:../../org/pyretechnics.org::suppression-imports][suppression-imports]]
import cython
import cython as cy
import numpy as np
if cython.compiled:
    from cython.cimports.numpy import ndarray
    from cython.cimports.libc.math import sqrt
    from cython.cimports.pyretechnics.cy_types import pyidx
else:
    from numpy import ndarray
    from math import sqrt
    from pyretechnics.py_types import pyidx
# suppression-imports ends here
# [[file:../../org/pyretechnics.org::order-frontier-cells][order-frontier-cells]]
def get_ordered_sides(y: cy.int, x: cy.int, is_burned: cy.bint) -> set:
    north: cy.float = y + 0.5
    south: cy.float = y - 0.5
    east : cy.float = x + 0.5
    west : cy.float = x - 0.5

    # Define cell sides as (y1,x1,y2,x2) where (y1,x1) and (y2,x2) are cell corners and (y1,x1)->(y2,x2)
    # is a directional cell side with an unburned cell to the left and a burned cell to the right.
    if is_burned:
        # Clockwise
        north_side: tuple[cy.float, cy.float, cy.float, cy.float] = (north, west, north, east)
        east_side : tuple[cy.float, cy.float, cy.float, cy.float] = (north, east, south, east)
        south_side: tuple[cy.float, cy.float, cy.float, cy.float] = (south, east, south, west)
        west_side : tuple[cy.float, cy.float, cy.float, cy.float] = (south, west, north, west)
    else:
        # Counterclockwise
        north_side: tuple[cy.float, cy.float, cy.float, cy.float] = (north, east, north, west)
        west_side : tuple[cy.float, cy.float, cy.float, cy.float] = (north, west, south, west)
        south_side: tuple[cy.float, cy.float, cy.float, cy.float] = (south, west, south, east)
        east_side : tuple[cy.float, cy.float, cy.float, cy.float] = (south, east, north, east)

    return {north_side, east_side, south_side, west_side}


def get_frontier_sides(frontier_cells: dict) -> set:
    burned_cells       : tuple[ndarray, ndarray] = frontier_cells["burned_cells"]
    unburned_cells     : tuple[ndarray, ndarray] = frontier_cells["unburned_cells"]
    burned_ys          : cy.int[::1]             = burned_cells[0]
    burned_xs          : cy.int[::1]             = burned_cells[1]
    unburned_ys        : cy.int[::1]             = unburned_cells[0]
    unburned_xs        : cy.int[::1]             = unburned_cells[1]
    burned_cell_sides  : set                     = set()
    unburned_cell_sides: set                     = set()
    i                  : pyidx
    for i in range(len(burned_ys)):
        burned_cell_sides.update(get_ordered_sides(burned_ys[i], burned_xs[i], True))

    for i in range(len(unburned_ys)):
        unburned_cell_sides.update(get_ordered_sides(unburned_ys[i], unburned_xs[i], False))

    frontier_sides: set = set.intersection(burned_cell_sides, unburned_cell_sides)

    return frontier_sides


def get_frontier_corner_graph(frontier_sides):
    frontier_corner_graph = {}

    for side in frontier_sides:
        (y1, x1, y2, x2) = side
        corner1          = (y1, x1)
        corner2          = (y2, x2)
        # NOTE: If corner1 has multiple corner2 values, we retain both.
        prev_corner2 = frontier_corner_graph.get(corner1)
        if prev_corner2:
            frontier_corner_graph[corner1] = [prev_corner2, corner2]
        else:
            frontier_corner_graph[corner1] = corner2

    return frontier_corner_graph


def group_corners_by_type(frontier_corner_graph):
    parent_corners = set(frontier_corner_graph.keys())
    child_corners  = set()
    for corner in frontier_corner_graph.values():
        if isinstance(corner, list):
            child_corners.update(corner)
        else:
            child_corners.add(corner)
    root_corners = set.difference(parent_corners, child_corners)
    leaf_corners = set.difference(child_corners, parent_corners)
    return {
        "parent_corners": parent_corners,
        "child_corners" : child_corners,
        "root_corners"  : root_corners,
        "leaf_corners"  : leaf_corners,
    }


def get_side_direction(corner1, corner2):
    (y1, x1) = corner1
    (y2, x2) = corner2
    if y1 == y2:
        if x2 > x1:
            return "east"
        else:
            return "west"
    else:
        if y2 > y1:
            return "north"
        else:
            return "south"


get_left_direction = {
    "north": "west",
    "east" : "north",
    "south": "east",
    "west" : "south",
}


def visited_corner_before(corner, visited_corners):
    if isinstance(corner, list):
        return False
    else:
        return corner in visited_corners


def iterative_graph_traversal(frontier_corner_graph, root_corner):
    corner_sequence    = [root_corner]
    visited_corners    = set(corner_sequence)
    grandparent_corner = None
    parent_corner      = root_corner
    child_corner       = frontier_corner_graph.get(parent_corner)
    while(child_corner and not visited_corner_before(child_corner, visited_corners)):
        if isinstance(child_corner, list):
            (child_corner1, child_corner2) = child_corner
            if grandparent_corner is None:
                # NOTE: Only happens if the root_corner branches toward two child_corners (should be rare).
                # Use child_corner1 in the absence of better information
                corner_sequence.append(child_corner1)
                visited_corners.add(child_corner1)
                grandparent_corner = parent_corner
                parent_corner      = child_corner1
                child_corner       = frontier_corner_graph.get(parent_corner)
            else:
                # NOTE: Prefer a left turn over a right turn as this should increase polygon convexity
                #       since we are traversing the frontier in a clockwise direction.
                incoming_direction      = get_side_direction(grandparent_corner, parent_corner)
                outgoing_direction1     = get_side_direction(parent_corner, child_corner1)
                outgoing_direction2     = get_side_direction(parent_corner, child_corner2)
                outgoing_left_direction = get_left_direction[incoming_direction]
                if outgoing_left_direction == outgoing_direction1:
                    # Use child_corner1
                    corner_sequence.append(child_corner1)
                    visited_corners.add(child_corner1)
                    grandparent_corner = parent_corner
                    parent_corner      = child_corner1
                    child_corner       = frontier_corner_graph.get(parent_corner)
                else:
                    # Use child_corner2
                    corner_sequence.append(child_corner2)
                    visited_corners.add(child_corner2)
                    grandparent_corner = parent_corner
                    parent_corner      = child_corner2
                    child_corner       = frontier_corner_graph.get(parent_corner)
        else:
            corner_sequence.append(child_corner)
            visited_corners.add(child_corner)
            grandparent_corner = parent_corner
            parent_corner      = child_corner
            child_corner       = frontier_corner_graph.get(parent_corner)
    if (child_corner and visited_corner_before(child_corner, visited_corners)):
        # Ended on a cycle, so add it to corner_sequence
        corner_sequence.append(child_corner)
    return corner_sequence


def corners_to_ordered_frontier_cells(corner_sequence: list) -> dict:
    num_sides  : cy.int      = len(corner_sequence) - 1
    burned_ys  : cy.int[::1] = np.zeros(num_sides, dtype=np.int32)
    burned_xs  : cy.int[::1] = np.zeros(num_sides, dtype=np.int32)
    unburned_ys: cy.int[::1] = np.zeros(num_sides, dtype=np.int32)
    unburned_xs: cy.int[::1] = np.zeros(num_sides, dtype=np.int32)
    i          : pyidx
    for i in range(num_sides):
        parent_corner: tuple[cy.float, cy.float] = corner_sequence[i]
        child_corner : tuple[cy.float, cy.float] = corner_sequence[i+1]
        y1           : cy.float                  = parent_corner[0]
        x1           : cy.float                  = parent_corner[1]
        y2           : cy.float                  = child_corner[0]
        x2           : cy.float                  = child_corner[1]
        direction    : str                       = get_side_direction(parent_corner, child_corner)
        y            : cy.int
        y_up         : cy.int
        y_down       : cy.int
        x            : cy.int
        x_left       : cy.int
        x_right      : cy.int
        if direction == "north":
            y              = int(y1 + 0.5)
            x_left         = int(x1 - 0.5)
            x_right        = int(x1 + 0.5)
            burned_ys[i]   = y             # east_cell
            burned_xs[i]   = x_right       # east_cell
            unburned_ys[i] = y             # west_cell
            unburned_xs[i] = x_left        # west_cell
        elif direction == "east":
            y_up           = int(y1 + 0.5)
            y_down         = int(y1 - 0.5)
            x              = int(x1 + 0.5)
            burned_ys[i]   = y_down        # south_cell
            burned_xs[i]   = x             # south_cell
            unburned_ys[i] = y_up          # north_cell
            unburned_xs[i] = x             # north_cell
        elif direction == "south":
            y              = int(y2 + 0.5)
            x_left         = int(x1 - 0.5)
            x_right        = int(x1 + 0.5)
            burned_ys[i]   = y             # west_cell
            burned_xs[i]   = x_left        # west_cell
            unburned_ys[i] = y             # east_cell
            unburned_xs[i] = x_right       # east_cell
        else: # direction == "west"
            y_up           = int(y1 + 0.5)
            y_down         = int(y1 - 0.5)
            x              = int(x2 + 0.5)
            burned_ys[i]   = y_up          # north_cell
            burned_xs[i]   = x             # north_cell
            unburned_ys[i] = y_down        # south_cell
            unburned_xs[i] = x             # south_cell
    return {
        "burned_cells"  : (burned_ys, burned_xs),
        "unburned_cells": (unburned_ys, unburned_xs),
    }


def order_frontier_cells(frontier_cells):
    frontier_sides         = get_frontier_sides(frontier_cells)
    frontier_corner_graph  = get_frontier_corner_graph(frontier_sides)
    corner_info            = group_corners_by_type(frontier_corner_graph)
    ordered_frontier_cells = []

    # Process all linestrings
    for root_corner in corner_info["root_corners"]:
        corner_sequence = iterative_graph_traversal(frontier_corner_graph, root_corner)
        for corner in corner_sequence:
            frontier_corner_graph.pop(corner, None)
        ordered_frontier_cells.append(corners_to_ordered_frontier_cells(corner_sequence))

    # Process all rings
    while(len(frontier_corner_graph) > 0):
        parent_corner   = next(iter(frontier_corner_graph))
        corner_sequence = iterative_graph_traversal(frontier_corner_graph, parent_corner)
        for corner in corner_sequence:
            frontier_corner_graph.pop(corner, None)
        ordered_frontier_cells.append(corners_to_ordered_frontier_cells(corner_sequence))

    return ordered_frontier_cells
# order-frontier-cells ends here
# [[file:../../org/pyretechnics.org::select-frontier-cells-for-suppression][select-frontier-cells-for-suppression]]
def get_frontier_distances(frontier, cell_height, cell_width):
    burned_cells           = frontier["burned_cells"]
    burned_ys              = burned_cells[0]
    burned_xs              = burned_cells[1]
    unburned_cells         = frontier["unburned_cells"]
    unburned_ys            = unburned_cells[0]
    unburned_xs            = unburned_cells[1]
    frontier_length        = len(burned_ys)
    distances              = np.zeros(frontier_length, dtype=np.float32)
    half_diagonal_distance = sqrt(cell_height * cell_height + cell_width * cell_width)/2.0

    # Calculate the frontier length between each cell pair
    for i in range(frontier_length):
        if (burned_ys[i] == unburned_ys[i]):
            # vertical cell side
            distances[i] = cell_height
        else:
            # horizonal cell side
            distances[i] = cell_width

    # If two cell pairs share the same unburned cell, we reduce their respective distances
    # such that their sum will equal the length of a diagonal line across the unburned cell
    for i in range(frontier_length - 1):
        if unburned_ys[i] == unburned_ys[i+1] and unburned_xs[i] == unburned_xs[i+1]:
            distances[i] = half_diagonal_distance

    for i in range(1, frontier_length):
        if unburned_ys[i] == unburned_ys[i-1] and unburned_xs[i] == unburned_xs[i-1]:
            distances[i] = half_diagonal_distance

    if unburned_ys[frontier_length-1] == unburned_ys[0] and unburned_xs[frontier_length-1] == unburned_xs[0]:
        distances[frontier_length-1] = half_diagonal_distance
        distances[0]                 = half_diagonal_distance

    return distances


def score_frontier_cells(ordered_frontier_cells, cell_height, cell_width,
                         suppression_priority_function, fireline_construction_rate_function):
    frontier_cell_scores   = []

    for frontier in ordered_frontier_cells:
        # Calculate the suppression priority of each cell pair
        scores = suppression_priority_function(frontier)

        # Calculate the fireline construction rate of each cell pair
        rates = fireline_construction_rate_function(frontier)

        # Calculate the fireline length required for each cell pair
        distances = get_frontier_distances(frontier, cell_height, cell_width)

        # Calculate the fireline construction times between each cell pair
        times = distances / rates

        frontier_cell_scores.append({
            "scores"   : scores,
            "rates"    : rates,
            "distances": distances,
            "times"    : times,
        })

    return frontier_cell_scores


def select_frontier_cells_for_suppression(ordered_frontier_cells, frontier_cell_scores, max_construction_time):
    selected_frontier_cells = ()
    top_path_score          = 0.0

    for i in range(len(ordered_frontier_cells)):
        frontier        = ordered_frontier_cells[i]
        unburned_cells  = frontier["unburned_cells"]
        unburned_ys     = unburned_cells[0]
        unburned_xs     = unburned_cells[1]
        frontier_length = len(unburned_ys)
        frontier_scores = frontier_cell_scores[i]
        scores          = frontier_scores["scores"]
        times           = frontier_scores["times"]

        # Use a moving window to find the sequence of frontier cells with the maximum total score
        # whose total frontier length can be contained within max_construction_time
        start_idx         = 0
        stop_idx          = 0
        path_score        = 0.0
        construction_time = 0.0
        best_start_idx    = 0
        best_stop_idx     = 0
        best_path_score   = 0.0

        while(stop_idx < frontier_length):
            # Expand the moving window until we reach max_construction_time or the end of the frontier
            while(stop_idx < frontier_length and construction_time + times[stop_idx] < max_construction_time):
                path_score        += scores[stop_idx]
                construction_time += times[stop_idx]
                stop_idx          += 1

            # Record this path's info if its score is the highest one that we've seen yet
            if path_score > best_path_score:
                best_start_idx  = start_idx   # inclusive
                best_stop_idx   = stop_idx    # exclusive
                best_path_score = path_score

            # Shift the moving window forward by one step
            path_score        -= scores[start_idx]
            construction_time -= times[start_idx]
            start_idx         += 1

            # Detect if stop_idx is not moving forward and jump over the cell that is blocking fireline growth
            if start_idx > stop_idx:
                stop_idx = start_idx

        if best_path_score > top_path_score:
            top_path_score          = best_path_score
            selected_ys             = unburned_ys[best_start_idx:best_stop_idx]
            selected_xs             = unburned_xs[best_start_idx:best_stop_idx]
            selected_frontier_cells = (selected_ys, selected_xs)

    return selected_frontier_cells
# select-frontier-cells-for-suppression ends here
# [[file:../../org/pyretechnics.org::build-firelines][build-firelines]]
def build_firelines(fuel_model_matrix, frontier_cells, cell_height, cell_width, suppression_priority_function,
                    fireline_construction_rate_function, max_construction_time, nonburnable_fuel_model=91):
    # Order the frontier cells into a list of independent, spatially contiguous perimeters
    ordered_frontier_cells = order_frontier_cells(frontier_cells)

    # Assign suppression priorities and fireline construction times to each perimeter's cells
    frontier_cell_scores = score_frontier_cells(ordered_frontier_cells,
                                                cell_height,
                                                cell_width,
                                                suppression_priority_function,
                                                fireline_construction_rate_function)

    # Select cells to suppress for each perimeter
    selected_frontier_cells = select_frontier_cells_for_suppression(ordered_frontier_cells,
                                                                    frontier_cell_scores,
                                                                    max_construction_time)

    # Change the values in fuel_model_matrix to nonburnable_fuel_model for all selected_frontier_cells
    fuel_model_matrix[selected_frontier_cells] = nonburnable_fuel_model

    return fuel_model_matrix
# build-firelines ends here
# [[file:../../org/pyretechnics.org::example-suppression-priority-functions][example-suppression-priority-functions]]
def make_sdi_suppression_priority_fn(sdi_matrix):
    max_sdi = np.max(sdi_matrix) + 1.0

    def sdi_suppression_priority(frontier):
        return max_sdi - sdi_matrix[frontier["unburned_cells"]]

    return sdi_suppression_priority


def make_flame_length_suppression_priority_fn(flame_length_matrix):
    max_flame_length = np.max(flame_length_matrix) + 1.0

    def flame_length_suppression_priority(frontier):
        return max_flame_length - flame_length_matrix[frontier["burned_cells"]]

    return flame_length_suppression_priority
# example-suppression-priority-functions ends here
