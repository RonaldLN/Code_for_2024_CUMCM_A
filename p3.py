from utils import *
import numpy as np

def calculate_spiral_positions_for_pitch(pitch, distance_from_center):
    """
    Calculate positions on the spiral given a pitch and a specific distance from the center.
    
    Args:
        pitch (float): The pitch value to calculate positions for.
        distance_from_center (float): Distance from the center to place the first node (front handle of the dragon head).
        
    Returns:
        list of tuples: The positions of the nodes on the spiral.
    """
    # Example spiral calculation logic (you need to replace this with the actual formula for your spiral)
    theta = 2 * np.pi * distance_from_center / pitch
    positions = calculate_spiral_positions(theta, pitch=pitch)
    return positions

def check_collision_for_pitch(positions, width, node_to_edge_distance):
    """
    Check if any two rectangles collide for a given set of positions.
    
    Args:
        positions (list of tuples): List of positions representing the spiral nodes.
        width (float): Width of the rectangles.
        node_to_edge_distance (float): Distance from node to rectangle edge.
        
    Returns:
        bool: True if a collision is detected, False otherwise.
    """
    # Define nodes for head rectangle (rect1) and first body rectangle (rect2)
    head_rect_nodes = [positions[0][:2], positions[1][:2]]
    first_body_rect_nodes = [positions[1][:2], positions[2][:2]]
    
    # Define range for further body parts checking
    theta_begin = positions[0][2] + 4 / 3 * np.pi
    theta_end = positions[1][2] + 8 / 3 * np.pi

    # Iterate through remaining nodes to check collisions
    for i in range(2, len(positions) - 1):
        if positions[i][2] < theta_begin:
            continue
        elif positions[i][2] > theta_end:
            break
        # Define nodes of the current body rectangle
        current_rect_nodes = [positions[i][:2], positions[i + 1][:2]]
        
        # Check for collision between head and body parts
        if check_collision(head_rect_nodes, current_rect_nodes, width, node_to_edge_distance) or \
           check_collision(first_body_rect_nodes, current_rect_nodes, width, node_to_edge_distance):
            return True  # Collision detected
    return False  # No collision detected

def find_collision_pitch(width, node_to_edge_distance, distance_from_center, initial_pitch=0.55, pitch_step=0.1, precision=1e-6):
    """
    Find the pitch value at which a collision occurs.
    
    Args:
        width (float): Width of the rectangles.
        node_to_edge_distance (float): Distance from node to rectangle edge.
        distance_from_center (float): Distance of the first node from the center of the spiral.
        initial_pitch (float): Starting pitch value to check for collisions.
        pitch_step (float): Step size to increment the pitch.
        precision (float): Desired precision for finding the exact collision pitch.
    
    Returns:
        float: The precise pitch at which a collision occurs.
    """
    pitch = initial_pitch
    step = pitch_step

    while step > precision:
        collision_occurred = False
        while not collision_occurred:
            # Calculate the positions of nodes at the current pitch
            positions = calculate_spiral_positions_for_pitch(pitch, distance_from_center)
            
            # Check if there is a collision at this pitch
            if check_collision_for_pitch(positions, width, node_to_edge_distance):
                collision_occurred = True
            else:
                # Increment the pitch by the step size
                pitch -= step
        
        # Reduce the step size for more precision and adjust the pitch
        pitch += step  # Go back to the start of the interval where collision was first detected
        step /= 10  # Reduce step size for finer checking

    return pitch

if __name__ == "__main__":
    # Example usage
    width = 0.3  # Width of the rectangles
    distance_from_center = 4.5  # Distance from the center for the first node
    node_to_edge_distance = 0.275  # Distance from node to rectangle edge
    
    # Find collision pitch starting from 0.55
    collision_pitch = find_collision_pitch(width=width, node_to_edge_distance=node_to_edge_distance, distance_from_center=distance_from_center)
    print(f"Collision occurs at approximately pitch = {collision_pitch:.6f} meters.")
    print(f"theta_s = {2 * np.pi * distance_from_center / collision_pitch} radians")
