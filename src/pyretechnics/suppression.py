# [[file:../../org/pyretechnics.org::order-frontier-cells][order-frontier-cells]]
def get_ordered_sides(cell_index, is_burned):
    y     = cell_index[0]
    x     = cell_index[1]
    north = y + 0.5
    south = y - 0.5
    east  = x + 0.5
    west  = x - 0.5

    # Define cell sides as (y1,x1,y2,x2) where (y1,x1) and (y2,x2) are cell corners and (y1,x1)->(y2,x2)
    # is a directional cell side with an unburned cell to the left and a burned cell to the right.
    if is_burned:
        # Clockwise
        north_side = (north, west, north, east)
        east_side  = (north, east, south, east)
        south_side = (south, east, south, west)
        west_side  = (south, west, north, west)
    else:
        # Counterclockwise
        north_side = (north, east, north, west)
        west_side  = (north, west, south, west)
        south_side = (south, west, south, east)
        east_side  = (south, east, north, east)

    return {north_side, east_side, south_side, west_side}


def get_frontier_sides(frontier_cells):
    burned_cells        = frontier_cells["burned_cells"]
    unburned_cells      = frontier_cells["unburned_cells"]
    burned_cell_sides   = set()
    unburned_cell_sides = set()

    for cell in burned_cells:
        burned_cell_sides.update(get_ordered_sides(cell, True))

    for cell in unburned_cells:
        unburned_cell_sides.update(get_ordered_sides(cell, False))

    frontier_sides = set.intersection(burned_cell_sides, unburned_cell_sides)

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


def iterative_graph_traversal(frontier_corner_graph, root_corner):
    corner_sequence    = [root_corner]
    visited_corners    = set(corner_sequence)
    grandparent_corner = None
    parent_corner      = root_corner
    child_corner       = frontier_corner_graph.get(parent_corner)
    while(child_corner and child_corner not in visited_corners):
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
                outgoing_left_direction = get_left_direction(incoming_direction)
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
    if (child_corner and child_corner in visited_corners):
        # Ended on a cycle, so add it to corner_sequence
        corner_sequence.append(child_corner)
    return corner_sequence


def corners_to_ordered_frontier_cells(corner_sequence):
    burned_cells   = []
    unburned_cells = []
    parent_corner  = corner_sequence[0]
    for child_corner in corner_sequence[1:]:
        (y1, x1)  = parent_corner
        (y2, x2)  = child_corner
        direction = get_side_direction(parent_corner, child_corner)
        if direction == "north":
            y         = int(y1 + 0.5)
            x_left    = int(x1 - 0.5)
            x_right   = int(x1 + 0.5)
            east_cell = (y, x_right)
            west_cell = (y, x_left)
            burned_cells.append(east_cell)
            unburned_cells.append(west_cell)
        elif direction == "east":
            y_up       = int(y1 + 0.5)
            y_down     = int(y1 - 0.5)
            x          = int(x1 + 0.5)
            north_cell = (y_up  , x)
            south_cell = (y_down, x)
            burned_cells.append(south_cell)
            unburned_cells.append(north_cell)
        elif direction == "south":
            y         = int(y2 + 0.5)
            x_left    = int(x1 - 0.5)
            x_right   = int(x1 + 0.5)
            east_cell = (y, x_right)
            west_cell = (y, x_left)
            burned_cells.append(west_cell)
            unburned_cells.append(east_cell)
        else: # direction == "west"
            y_up       = int(y1 + 0.5)
            y_down     = int(y1 - 0.5)
            x          = int(x2 + 0.5)
            north_cell = (y_up  , x)
            south_cell = (y_down, x)
            burned_cells.append(north_cell)
            unburned_cells.append(south_cell)
        parent_corner = child_corner
    return {
        "burned_cells"  : burned_cells,
        "unburned_cells": unburned_cells,
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
# [[file:../../org/pyretechnics.org::#selecting-frontier-cells-for-suppression][Selecting Frontier Cells for Suppression:1]]
import numpy as np


def score_frontier_cells(ordered_frontier_cells, scoring_function):
    frontier_cell_scores = []

    for frontier in ordered_frontier_cells:
        burned_cells           = frontier["burned_cells"]
        unburned_cells         = frontier["unburned_cells"]
        frontier_length        = len(burned_cells)
        scores                 = np.zeros(frontier_length)
        distances              = np.ones(frontier_length)
        half_diagonal_distance = 0.7071 # approximately math.sqrt(2)/2

        # Score each cell pair
        for i in range(frontier_length):
            scores[i] = scoring_function(burned_cells[i], unburned_cells[i])

        # Calculate the distances along the frontier
        # NOTE: If two cell pairs share the same unburned cell, we reduce their respective distances
        #       such that their sum will equal the length of a diagonal line across the unburned cell.
        for i in range(frontier_length - 1):
            if unburned_cells[i] == unburned_cells[i+1]:
                distances[i] = half_diagonal_distance

        for i in range(1, frontier_length):
            if unburned_cells[i] == unburned_cells[i-1]:
                distances[i] = half_diagonal_distance

        if unburned_cells[frontier_length-1] == unburned_cells[0]:
            distances[frontier_length-1] = half_diagonal_distance
            distances[0]                 = half_diagonal_distance

        frontier_cell_scores.append({
            "scores"   : scores,
            "distances": distances,
        })

    return frontier_cell_scores


# NOTE: max_suppression_length should be expressed in cell side lengths
def select_frontier_cells_for_suppression(ordered_frontier_cells, frontier_cell_scores, max_suppression_length):
    if max_suppression_length < 1.0:
        raise ValueError("The max_suppression_length must be at least 1.0 in order to be able "
                         + " to suppress at least one frontier cell.")

    selected_frontier_cells = []

    for i in range(len(ordered_frontier_cells)):
        frontier        = ordered_frontier_cells[i]
        unburned_cells  = frontier["unburned_cells"]
        frontier_length = len(unburned_cells)
        frontier_scores = frontier_cell_scores[i]
        scores          = frontier_scores["scores"]
        distances       = frontier_scores["distances"]

        # Use a moving window to find the sequence of frontier cells with the maximum total score
        # whose total frontier length is less than or equal to max_suppression_length.
        start_idx          = 0
        stop_idx           = 0
        path_score         = 0.0
        suppression_length = 0.0
        best_start_idx     = 0
        best_stop_idx      = 0
        best_path_score    = 0.0

        while(stop_idx < frontier_length):
            # Expand the moving window until we reach max_suppression_length or the end of the frontier
            while(stop_idx < frontier_length and suppression_length + distances[stop_idx] < max_suppression_length):
                path_score         += scores[stop_idx]
                suppression_length += distances[stop_idx]
                stop_idx           += 1

            # Record this path's info if its score is the highest one that we've seen yet
            if path_score > best_path_score:
                best_start_idx  = start_idx   # inclusive
                best_stop_idx   = stop_idx    # exclusive
                best_path_score = path_score

            # Shift the moving window forward by one step
            path_score         -= scores[start_idx]
            suppression_length -= distances[start_idx]
            start_idx          += 1

        selected_frontier_cells.append(unburned_cells[best_start_idx:best_stop_idx])

    return selected_frontier_cells
# Selecting Frontier Cells for Suppression:1 ends here
