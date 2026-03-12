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


# TODO: Remove debugging print statements
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

    if len(root_corners) == 0:
        if len(leaf_corners) == 0:
            # Begin anywhere and traverse one or more cycles
            print("Only cycles detected.")
        else:
            # Begin anywhere and traverse one or more cycles with one or more branching nodes that terminate on leaves
            print("No roots but some trailing leaves exist.")
    else:
        if len(leaf_corners) == 0:
            # Begin on roots and traverse one or more cycles with branching nodes that continue the cycles
            print("No leaves but some trailing roots exist.")
        else:
            # Begin on roots and traverse one or more cycles with one or more branching nodes that terminate on leaves
            print("Both roots and leaves detected.")

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
    corner_sequence        = iterative_graph_traversal(frontier_corner_graph, list(corner_info["parent_corners"])[0]) # FIXME: Explore all paths from the root_corners until all corners have been visited. Return a list of lists.
    ordered_frontier_cells = corners_to_ordered_frontier_cells(corner_sequence)

    return ordered_frontier_cells
# order-frontier-cells ends here
