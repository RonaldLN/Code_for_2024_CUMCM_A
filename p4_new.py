import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

from utils.p4_turn_path import *

PITCH = 1.7

def spiral_arc_length(theta_a, theta_b, pitch=PITCH):
    """
    Calculate the arc length of the spiral from theta_a to theta_b.
    """
    integrand = lambda theta: np.sqrt((pitch / (2 * np.pi))**2 + (pitch / (2 * np.pi) * theta)**2)
    arc_length, _ = quad(integrand, theta_a, theta_b)
    return arc_length

def solve_spiral_theta_for_s(s, theta_lower, pitch=PITCH):
    """
    Find the value of theta that corresponds to a specific arc length (distance) s.
    """
    equation = lambda theta: spiral_arc_length(theta_lower, theta, pitch) - s
    theta_solution = fsolve(equation, s)[0]
    return theta_solution

def s_to_position(s, theta_start, theta_end, center1, center2, pitch, angle1, angle2, length1, length2):
    """
    Maps the parameter 's' to a position along the path.
    When s = 0, it returns the point corresponding to 'theta_start'.
    
    Parameters:
    - s: The parameter value, representing a distance along the path.
    - theta_start: Starting angle for the initial spiral.
    - theta_end: Ending angle for the initial spiral.
    - center1: Center of the first circular arc.
    - center2: Center of the second circular arc.
    - pitch: Pitch for the spiral.
    - angle1: Starting angle for the first circular arc.
    - angle2: Starting angle for the second circular arc.
    - length1: Length of the path segment for the first semi-circle.
    - length2: Length of the path segment for the second semi-circle.
    
    Returns:
    - position: The position (x, y) corresponding to the parameter 's'.
    """
    
    # # Determine total length for normalization and handling path transitions
    # spiral_length = np.abs(theta_end - theta_start)  # Approximation of spiral length
    # total_length = spiral_length + length1 + length2 + spiral_length

    point_a = theta_to_spiral_position(theta_start, pitch)
    point_b = rotate_around_center(point_a, center2, -angle2)
    point_c = theta_to_spiral_position(theta_end, -pitch)

    r1 = length1 / angle1
    r2 = 2 * r1

    # Handle s = 0 case
    if s <= 0:
        theta = solve_spiral_theta_for_s(-s, theta_start, pitch)
        return theta_to_spiral_position(theta, pitch)

    # Determine which path segment 's' falls into
    # if s <= spiral_length:
    #     # First spiral segment
    #     theta = theta_start + (s / spiral_length) * (theta_end - theta_start)
    #     return theta_to_spiral_position(theta, -pitch)
    if s <= length2:
        # First circular arc around center2
        s_arc = s
        angle = s_arc / r2
        return rotate_around_center(point_a, center2, -angle)
    elif s <= length1 + length2:
        # Second circular arc around center1
        s_arc = s - length2
        angle = s_arc / r1
        return rotate_around_center(point_c, center1, angle-angle1)
    else:
        # Final spiral segment
        s_spiral = s - length1 - length2
        theta = solve_spiral_theta_for_s(s_spiral, theta_end, -pitch)
        return theta_to_spiral_position(theta, -pitch)

# Helper function to compute position for a spiral based on theta
def theta_to_spiral_position(theta, pitch):
    """
    Calculate the (x, y) position on a spiral given an angle (theta) and pitch.
    
    Parameters:
    - theta: Angle in radians.
    - pitch: Pitch of the spiral.
    
    Returns:
    - (x, y): Coordinates of the point on the spiral.
    """
    r = pitch * theta  # Simple Archimedean spiral
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return (x, y)

# Helper function to rotate a point around a given center
def rotate_around_center(point, center, angle):
    """
    Rotate a point around a center by a given angle.
    
    Parameters:
    - point: (x, y) coordinates of the point.
    - center: (x, y) coordinates of the center.
    - angle: Angle in radians to rotate.
    
    Returns:
    - (x_rot, y_rot): Rotated (x, y) coordinates.
    """
    # def rotate_vector(vector, angle):
    #     if not isinstance(angle, np.ndarray):
    #         angle = np.array([angle])
    #     rotation_matrix = np.array([[np.cos(angle)[0], -np.sin(angle)[0]], [np.sin(angle)[0], np.cos(angle)[0]]])
    #     return np.dot(rotation_matrix, vector)

    # def rotate_point_on_circle(center, point, angle):
    #     relative_point = np.array(point) - np.array(center)
    #     rotated_point = rotate_vector(relative_point, angle)
    #     return np.array(rotated_point + np.array(center))
    x, y = point
    cx, cy = center
    x -= cx
    y -= cy
    x_rot = x * np.cos(angle) - y * np.sin(angle) + cx
    y_rot = x * np.sin(angle) + y * np.cos(angle) + cy
    return (x_rot, y_rot)
    # return rotate_point_on_circle(center, point, angle)

def distance_between_points(p1, p2):
    return np.linalg.norm(np.array(p2) - np.array(p1))

def arc_length(r, angle):
    return r * angle

NUM_NODES = 223
D1 = 2.86
D = 1.65

def calculate_positions_all_paths(t, theta_start, theta_end, center1, center2, angle1, angle2, pitch=PITCH):
    r1 = distance_between_points(center1, center2) / 3
    r2 = 2 * r1
    len1 = arc_length(r1, angle1)
    len2 = arc_length(r2, angle2)

    s = t

    def next_position(node_i, s_i, distance):
        def equation(s):
            return distance_between_points(node_i, s_to_position(s, theta_start, theta_end, center1, center2, pitch, angle1, angle2, len1, len2)) - distance
        
        s_guess = s_i - 4 * distance
        s_next = fsolve(equation, s_guess)[0]
        return s_next
    
    # Initialize the first node's position
    s_1 = s
    pos_1 = s_to_position(s_1, theta_start, theta_end, center1, center2, pitch, angle1, angle2, len1, len2)
    positions = [pos_1]

    # Calculate positions for subsequent nodes
    s_current = s_1
    node_current = pos_1

    # Loop through the nodes
    for i in range(NUM_NODES):
        # Use d1 for the first distance and d for all subsequent distances
        distance = D1 if i == 0 else D
        
        # Calculate the next s based on the current s and position
        s_next = next_position(node_current, s_current, distance)
        pos_next = s_to_position(s_next, theta_start, theta_end, center1, center2, pitch, angle1, angle2, len1, len2)
        
        # Append the new position
        positions.append(pos_next)
        
        # Update current s and position for the next iteration
        s_current = s_next
        node_current = pos_next

    return positions


if __name__ == "__main__":
    theta_start = 12
    optimal_theta_end = find_theta_end_for_r1(theta_start, pitch=PITCH)
    _, (_, angle1, angle2, _, center1, center2) = u_turn_path_length((theta_start, optimal_theta_end), pitch=PITCH)

    # positions = calculate_positions_all_paths(10, theta_start, optimal_theta_end, center1, center2, angle1, angle2, pitch=PITCH)
    # np.savetxt("positions.csv", positions, delimiter=",")

    r1 = distance_between_points(center1, center2) / 3

    positions = []
    for i in range(-1000, 1000):
        positions.append(s_to_position(i / 10, theta_start, optimal_theta_end, center1, center2, PITCH, angle1, angle2, arc_length(r1, angle1), arc_length(2 * r1, angle2)))
    np.savetxt("positions.csv", positions, delimiter=",")