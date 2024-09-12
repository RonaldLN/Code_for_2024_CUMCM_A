from .p1_velocity import calculate_tangent_direction
from .p4_position import theta_to_spiral_position
import numpy as np
from scipy.optimize import minimize, fsolve

def calculate_arc_center(position, tangent, radius):
    """
    Calculate the center of the arc that is perpendicular to the tangent at a given position.

    Args:
        position (tuple): The (x, y) position on the spiral.
        tangent (tuple): The tangent direction (dx, dy).
        radius (float): The radius of the arc.

    Returns:
        (float, float): The (x, y) coordinates of the arc center.
    """
    # Normalize the tangent vector
    tangent_norm = np.array(tangent) / np.linalg.norm(tangent)
    # Perpendicular vector to the tangent (rotated by 90 degrees)
    perp_vector = np.array([-tangent_norm[1], tangent_norm[0]])
    # Arc center is offset by the radius along the perpendicular direction
    arc_center = np.array(position) + perp_vector * radius
    return arc_center

def arc_length(radius, angle):
    """
    Calculate the length of an arc given its radius and angle.

    Args:
        radius (float): The radius of the arc.
        angle (float): The central angle of the arc in radians.

    Returns:
        float: The length of the arc.
    """
    return radius * angle

def u_turn_path_length(params, pitch):
    """
    Calculate the total U-turn path length for a given pair of theta_start and theta_end.

    Args:
        params (tuple): The (theta_start, theta_end) values.
        pitch (float): The pitch of the spiral.
        d1 (float): Initial distance from the center of the spiral.

    Returns:
        float: The total length of the U-turn path.
    """
    theta_start, theta_end = params
    
    # Calculate positions and tangents at the start and end points
    pos_start = theta_to_spiral_position(theta_start, pitch=pitch)
    pos_end = theta_to_spiral_position(theta_end, pitch=-pitch)
    
    tangent_start = calculate_tangent_direction(theta_start, pitch=pitch)
    tangent_end = calculate_tangent_direction(theta_end, pitch=-pitch)

    perp_vector_start = np.array([-tangent_start[1], tangent_start[0]])
    perp_vector_end = np.array([-tangent_end[1], tangent_end[0]])

    def equation(r1):
        center1 = pos_end + r1 * perp_vector_end
        center2 = pos_start + 2 * r1 * perp_vector_start
        return np.linalg.norm(center2 - center1) - 3 * r1
    
    # Assume r1 is the smaller radius and r2 is twice r1
    r1 = np.linalg.norm(np.array(pos_end) - np.array(pos_start)) / 5  # Approximate starting value

    # Solve for r1
    r1 = fsolve(equation, r1)[0]
    r2 = 2 * r1
    
    # Calculate arc centers
    center1 = calculate_arc_center(pos_end, tangent_end, r1)
    center2 = calculate_arc_center(pos_start, tangent_start, r2)
    
    # Calculate angles for the arcs
    angle1 = np.arccos(np.dot((np.array(center2) - np.array(center1)), (np.array(pos_end) - np.array(center1))) / (r1 * r1 * 3))
    angle2 = np.arccos(np.dot((np.array(center1) - np.array(center2)), (np.array(pos_start) - np.array(center2))) / (r1 * 3 * r2))
    
    if np.dot((center2 - center1) + (pos_end - center1), tangent_end) > 0:
        angle1 = 2 * np.pi - angle1
    if np.dot((center1 - center2) + (pos_start - center2), tangent_start) > 0:
        angle2 = 2 * np.pi - angle2

    # print(center1, center2)
    # print(np.linalg.norm(center2 - center1), 3 * r1)

    # Calculate arc lengths
    length1 = arc_length(r1, angle1)
    length2 = arc_length(r2, angle2)
    
    # Total path length is the sum of both arc lengths
    return length1 + length2, (r1, angle1, angle2, length1 + length2, center1, center2)

def find_theta_end_for_r1(theta_start, pitch, target_r1=0.725, initial_step=np.pi/10, precision=1e-6):
    """
    Find the theta_end for a given theta_start that results in r1 being close to target_r1.
    
    Args:
        theta_start (float): The starting theta.
        pitch (float): The pitch of the spiral.
        target_r1 (float): The target value for r1.
        initial_step (float): The initial step size for theta_end.
        precision (float): The acceptable difference from the target_r1.
        
    Returns:
        float: The optimal theta_end value.
    """
    # Initialize theta_end and step size
    theta_end = theta_start - np.pi  # Start with a reasonable guess
    step = initial_step

    while step > precision:
        optimal_theta_end_found = False

        while not optimal_theta_end_found:
            # Calculate r1 for the current theta_end
            _, (r1, _, _, _, _, _) = u_turn_path_length((theta_start, theta_end), pitch)
            
            if r1 > target_r1:
                optimal_theta_end_found = True
            else:
                # Increment theta_end by the step size
                theta_end += step
        
        # Reduce the step size for more precision and adjust theta_end
        theta_end -= step  # Go back to the start of the interval where the optimal theta_end was first detected
        step /= 10  # Reduce step size for finer checking

    return theta_end

if __name__ == "__main__":
    # Example usage
    width = 0.3  # Width of the rectangles
    node_to_edge_distance = 0.275  # Distance from node to rectangle edge
    pitch = 1.7  # Pitch of the spiral
    d1 = 2.2  # Initial distance from the center

    # Find the optimal starting and ending points for the U-turn
    optimal_thetas = find_optimal_u_turn_path(pitch)
    print(f"Optimal theta_start = {optimal_thetas[0]:.6f}, theta_end = {optimal_thetas[1]:.6f}")