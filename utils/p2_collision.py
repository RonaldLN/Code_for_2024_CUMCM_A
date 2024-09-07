from .p1_position import calculate_spiral_positions
from .p1_get_theta_at_time import find_theta_at_time
import numpy as np

def calculate_vertices(node1, node2, width, node_to_edge_distance):
    """
    Calculate the four vertices of a rectangle given two endpoints (nodes), width, 
    and the distance from the nodes to the rectangle edges.
    """
    # Calculate the direction vector between the two nodes
    dx, dy = node2[0] - node1[0], node2[1] - node1[1]
    length = np.sqrt(dx**2 + dy**2)
    
    # Calculate the perpendicular vector to the line segment between nodes
    px, py = -dy * (width / (2 * length)), dx * (width / (2 * length))
    
    # Calculate the unit direction vector
    ux, uy = dx / length, dy / length
    
    # Adjust nodes to the actual rectangle edges
    adjusted_node1 = (node1[0] - ux * node_to_edge_distance, node1[1] - uy * node_to_edge_distance)
    adjusted_node2 = (node2[0] + ux * node_to_edge_distance, node2[1] + uy * node_to_edge_distance)
    
    # Compute the four vertices of the rectangle
    vertex1 = (adjusted_node1[0] + px, adjusted_node1[1] + py)
    vertex2 = (adjusted_node1[0] - px, adjusted_node1[1] - py)
    vertex3 = (adjusted_node2[0] + px, adjusted_node2[1] + py)
    vertex4 = (adjusted_node2[0] - px, adjusted_node2[1] - py)
    
    return [vertex1, vertex2, vertex3, vertex4]

def is_point_within_rectangle(point, node1, node2, width):
    """
    Check if a point is within a rectangle defined by two nodes and width.
    """
    dx, dy = node2[0] - node1[0], node2[1] - node1[1]
    length = np.sqrt(dx**2 + dy**2)

    # Calculate the perpendicular distance from the point to the line segment (node1-node2)
    distance = abs(dy * point[0] - dx * point[1] + node2[0] * node1[1] - node2[1] * node1[0]) / length
    
    # Check if distance is within the half-width of the rectangle
    return distance <= (width / 2)

def check_collision(inner_rect_nodes, outer_rect_nodes, width, node_to_edge_distance):
    """
    Check if any of the vertices of inner_rect are within outer_rect.
    """
    inner_rect_vertices = calculate_vertices(inner_rect_nodes[0], inner_rect_nodes[1], width, node_to_edge_distance)
    
    # Check if any of the vertices of inner_rect are within outer_rect
    for vertex in inner_rect_vertices:
        if is_point_within_rectangle(vertex, outer_rect_nodes[0], outer_rect_nodes[1], width):
            return True
    return False

def find_collision_time(width, node_to_edge_distance, current_time=0.0, initial_step=1.0, precision=1e-6):
    """
    Find the collision time with an initial 1-second interval refinement.
    
    Args:
        width (float): Width of the rectangles.
        node_to_edge_distance (float): Distance from the node to the rectangle edge.
        current_time (float): Starting time to check for collisions.
        initial_step (float): Initial time step for coarse checking.
        precision (float): Desired precision for finding the exact collision time.
    
    Returns:
        float: The precise time at which a collision occurs.
    """
    step = initial_step

    while step > precision:
        collision_occurred = False
        while not collision_occurred:
            # Get the positions of the nodes at the current time
            theta_s = find_theta_at_time(current_time)
            
            # Calculate positions of dragon head (rect1) and first body (rect2)
            positions = calculate_spiral_positions(theta_s)
            
            # Nodes of head rectangle (rect1) and first body rectangle (rect2)
            head_rect_nodes = [positions[0][:2], positions[1][:2]]
            first_body_rect_nodes = [positions[1][:2], positions[2][:2]]

            # Define range for further body parts checking
            theta_begin = theta_s + 4 / 3 * np.pi
            theta_end = positions[1][2] + 8 / 3 * np.pi
            
            # Iterate through all body parts to check for collision
            it = iter(positions)
            theta = next(it)[2]
            node1, node2 = None, None
            while theta < theta_end:
                if theta < theta_begin:
                    theta = next(it)[2]
                    continue
                pos = next(it)
                node1, node2, theta = pos[:2], node1, pos[2]
                if node2 is None:
                    continue
                elif check_collision(head_rect_nodes, [node1, node2], width, node_to_edge_distance) \
                    or check_collision(first_body_rect_nodes, [node1, node2], width, node_to_edge_distance):
                    collision_occurred = True
                    break

            if collision_occurred:
                break
            # Increment the time by step
            current_time += step
        
        # Reduce the step size for more precision and adjust the current time
        current_time -= step  # Go back to the start of the interval where collision first detected
        step /= 10  # Reduce step size for finer checking

    return current_time

if __name__ == "__main__":
    # Example usage
    width = 0.3  # Assume the width of the rectangles
    collision_time = find_collision_time(width=width, node_to_edge_distance=0.275)
    print(f"Collision occurs at approximately {collision_time:.6f} seconds.")
    
    # Start from a different time point with finer step
    collision_time = find_collision_time(width=width, node_to_edge_distance=0.275, current_time=300.0, initial_step=0.1)
    print(f"Collision occurs at approximately {collision_time:.6f} seconds.")
