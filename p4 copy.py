from utils import *
import numpy as np
from scipy.optimize import fsolve

from p4_new import calculate_positions_all_paths

PITCH = 1.7
NUM_NODES = 223
D1 = 2.86
D = 1.65
WIDTH = 0.3
NODE_TO_EDGE_DISTANCE = 0.275

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

def arc_length(r, angle):
    """
    Calculate the arc length of an arc with radius r and angle in radians.
    """
    return r * angle

def distance_between_points(p1, p2):
    """
    Calculate the Euclidean distance between two points.
    """
    return np.linalg.norm(np.array(p2) - np.array(p1))

def rotate_vector(vector, angle):
    """
    Rotate a 2D vector by a given angle in radians.
    """
    if not isinstance(angle, np.ndarray):
        angle = np.array([angle])
    rotation_matrix = np.array([[np.cos(angle)[0], -np.sin(angle)[0]], [np.sin(angle)[0], np.cos(angle)[0]]])
    return np.dot(rotation_matrix, vector)

def rotate_point_on_circle(center, point, angle, direction=1):
    """
    Rotate a point around a center by a given angle in radians.
    """
    if direction * angle < 0:
        return center
    relative_point = np.array(point) - np.array(center)
    rotated_point = rotate_vector(relative_point, angle)
    return np.array(rotated_point + np.array(center))

def calculate_spiral_theta_from_point(point, pitch=PITCH):
    """
    Calculate the angle theta of a point on the spiral.
    """
    r = np.linalg.norm(np.array(point))
    return 2 * np.pi * r / pitch

# def calculate_positions_all_paths(t, theta_start, theta_end, center1, center2, angle1, angle2, pitch=PITCH):
#     # if t < 0:
#     #     theta_s = solve_spiral_theta_for_s(-t, theta_start, pitch)
#     #     return calculate_spiral_positions(theta_s, pitch=PITCH)
    
#     # if t < arc_length(r2, angle2):
    
#     r1 = distance_between_points(center1, center2) / 3
#     r2 = 2 * r1
#     len1 = arc_length(r1, angle1)
#     len2 = arc_length(r2, angle2)

#     positions = []
#     node_index = 0
    
#     if t > len2 + len1:
#         theta_1 = solve_spiral_theta_for_s(t - len2 - len1, theta_end, -pitch)
#         node_1 = theta_to_spiral_position(theta_1, pitch=-pitch)
#         positions.append(node_1)
#         calculate_positions_stage_4(node_1, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
#     elif t > len2:
#         angle = (t - len2) / r1
#         node_1 = rotate_point_on_circle(center2, center1, angle)
#         positions.append(node_1)
#         calculate_positions_stage_3(node_1, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
#     elif t > 0:
#         angle = t / r2
#         node_1 = rotate_point_on_circle(center1, center2, -angle)
#         positions.append(node_1)
#         calculate_positions_stage_2(node_1, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
#     elif t == 0:
#         positions.append(theta_to_spiral_position(theta_start, pitch))
#         calculate_positions_stage_1(theta_to_spiral_position(theta_start, pitch), node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
#     else:
#         theta_1 = solve_spiral_theta_for_s(-t, theta_start, pitch)
#         node_1 = theta_to_spiral_position(theta_1, pitch)
#         positions.append(node_1)
#         calculate_positions_stage_1(node_1, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
    
#     return positions

def calculate_positions_stage_4(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch=PITCH):
    if node_index == NUM_NODES - 1:
        return positions
    
    d = D1 if node_index == 0 else D
    bound_point = theta_to_spiral_position(theta_end, -pitch)
    if distance_between_points(last_node, bound_point) < d:
        return calculate_positions_stage_3(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)

    def equation(theta):
        node = theta_to_spiral_position(theta, -pitch)
        return distance_between_points(last_node, node) - d
    
    theta_solution = fsolve(equation, theta_end)[0]
    node = theta_to_spiral_position(theta_solution, -pitch)
    positions.append(node)
    return calculate_positions_stage_4(node, node_index + 1, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)

def calculate_positions_stage_3(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch=PITCH):
    if node_index == NUM_NODES - 1:
        return positions
    
    d = D1 if node_index == 0 else D
    bound_point = (np.array(center2) - np.array(center1)) / 3 + np.array(center1)

    def equation(angle):
        node = rotate_point_on_circle(center1, bound_point, angle)
        return distance_between_points(last_node, node) - d
    
    solution = fsolve(equation, 0)
    if len(solution) == 0 or solution[0] <= 0 or solution[0] > angle1:
        return calculate_positions_stage_2(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
    angle_solution = solution[0]
    node = rotate_point_on_circle(center1, bound_point, angle_solution)
    positions.append(node)
    return calculate_positions_stage_3(node, node_index + 1, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)

def calculate_positions_stage_2(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch=PITCH):
    if node_index == NUM_NODES - 1:
        return positions
    
    d = D1 if node_index == 0 else D
    # bound_point = (np.array(center1) - np.array(center2)) * 2 / 3 + np.array(center2)
    bound_point = theta_to_spiral_position(theta_start, pitch)

    def equation(angle):
        node = rotate_point_on_circle(center2, bound_point, angle, direction=-1)
        return distance_between_points(last_node, node) - d
    
    solution = fsolve(equation, 0)
    if len(solution) == 0 or -solution[0] <= 0 or -solution[0] > angle2:
        return calculate_positions_stage_1(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)
    angle_solution = solution[0]
    node = rotate_point_on_circle(center2, bound_point, angle_solution)
    positions.append(node)
    return calculate_positions_stage_2(node, node_index + 1, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)

def calculate_positions_stage_1(last_node, node_index, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch=PITCH):
    if node_index == NUM_NODES - 1:
        return positions
    
    d = D1 if node_index == 0 else D
    def equation(theta):
        node = theta_to_spiral_position(theta, pitch)
        return distance_between_points(last_node, node) - d
    
    last_theta = calculate_spiral_theta_from_point(last_node, pitch)
    last_r = np.linalg.norm(np.array(last_node))
    theta_solution = fsolve(equation, max(theta_start, last_theta + d / last_r))[0]
    node = theta_to_spiral_position(theta_solution, pitch)
    positions.append(node)
    return calculate_positions_stage_1(node, node_index + 1, positions, theta_start, theta_end, center1, center2, angle1, angle2, pitch)

def check_collision_all_paths(positions, width=WIDTH, node_to_edge_distance=NODE_TO_EDGE_DISTANCE):
    head_rect_nodes = positions[:2]
    for i in range(2, len(positions) - 2):
        if check_collision(positions[i:i+2], head_rect_nodes, width, node_to_edge_distance):
            print(f"Collision detected at index {i}")  # Debugging output
            return True
    return False

def find_min_theta_start_no_collision(initial_theta_start, pitch=PITCH, theta_step=np.pi / 10, precision=1e-6):
    theta_start = initial_theta_start
    step = theta_step

    while step > precision:
        not_collision_occurred = False
        while not not_collision_occurred:
            optimal_theta_end = find_theta_end_for_r1(theta_start, pitch)
            _, (_, angle1, angle2, _, center1, center2) = u_turn_path_length((theta_start, optimal_theta_end), pitch)
            for t in range(0, 101):
                positions = calculate_positions_all_paths(t, theta_start, optimal_theta_end, center1, center2, angle1, angle2, pitch)
                if not check_collision_all_paths(positions):
                    not_collision_occurred = True
                    break
            
            if not_collision_occurred:
                break

            theta_start += step
        
        theta_start -= step
        step /= 10

    return theta_start

    
if __name__ == "__main__":
    theta_start = 12
    optimal_theta_end = find_theta_end_for_r1(theta_start, pitch=PITCH)
    _, (_, angle1, angle2, _, center1, center2) = u_turn_path_length((theta_start, optimal_theta_end), pitch=PITCH)
    # for t in range(-100, 101):
    #     # Calculate the positions of all nodes at time t
    #     positions = calculate_positions_all_paths(t, theta_start, optimal_theta_end, center1, center2, angle1, angle2, pitch=PITCH)
    #     print(f"Time: {t}s, Positions[0]: {positions[0]}")

    theta_start = find_min_theta_start_no_collision(8)
    print(f"Optimal theta_start: {theta_start}")

    # positions = calculate_positions_all_paths(56, theta_start, optimal_theta_end, center1, center2, angle1, angle2, pitch=PITCH)
    # np.savetxt("positions.csv", positions, delimiter=",")
