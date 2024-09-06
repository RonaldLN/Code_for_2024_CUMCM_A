import numpy as np
from .p1_position import calculate_spiral_positions  # Assuming p1_position.py contains the function to calculate positions

def calculate_tangent_direction(theta, pitch=0.55):
    """
    Calculate the unit tangent direction vector for a point on the Archimedean spiral.
    
    Args:
    - theta (float): Angle in radians.
    - pitch (float): Spiral pitch in meters per turn.
    
    Returns:
    - (vx, vy) (tuple of floats): Unit tangent vector (vx, vy) at angle theta.
    """
    dr_dtheta = pitch / (2 * np.pi)  # Derivative of r with respect to theta
    r = pitch / (2 * np.pi) * theta  # Radial distance for the given theta
    
    # Tangent vector components in Cartesian coordinates
    vx = dr_dtheta * np.cos(theta) - r * np.sin(theta)
    vy = dr_dtheta * np.sin(theta) + r * np.cos(theta)
    
    # Normalize the tangent vector to get a unit vector
    magnitude = np.sqrt(vx**2 + vy**2)
    vx /= magnitude
    vy /= magnitude
    
    return vx, vy

def calculate_velocities(positions, pitch=0.55, v1=1.0):
    """
    Calculate the velocities of all nodes based on the velocity of the first node.
    
    Args:
    - positions (list of tuples): List of (x, y, theta) positions of the nodes. 
                                  Each tuple is (x, y, theta) where theta is the angle in radians.
    - pitch (float): Spiral pitch in meters per turn.
    - v1 (float): Velocity of the first node (in m/s).
    
    Returns:
    - velocities (list of tuples): List of (vx, vy) velocity components of each node.
    """
    num_nodes = len(positions)
    velocities = [(0, 0)] * num_nodes  # Initialize velocity list with zero vectors
    
    # Calculate velocity for the first node (v1 = 1 m/s)
    theta_1 = positions[0][2]  # Extract theta for the first node
    vx_1, vy_1 = calculate_tangent_direction(theta_1, pitch)  # Get tangent direction for the first node
    velocities[0] = (v1 * vx_1, v1 * vy_1)  # Scale by the velocity v1
    
    # Calculate velocities for subsequent nodes
    for i in range(1, num_nodes):
        # Calculate direction of the line connecting node i to node i-1
        dx = positions[i][0] - positions[i-1][0]
        dy = positions[i][1] - positions[i-1][1]
        distance = np.sqrt(dx**2 + dy**2)  # Euclidean distance between nodes
        
        # Project the velocity of node i-1 along the line connecting to node i
        v_proj = (velocities[i-1][0] * dx + velocities[i-1][1] * dy) / distance

        # Tangent direction for node i
        theta_i = positions[i][2]  # Extract theta for node i
        vx_i, vy_i = calculate_tangent_direction(theta_i, pitch)  # Get tangent direction for node i
        
        # Calculate the velocity component in the direction of the tangent
        # Ensure the velocity maintains the component along the connecting line
        v_tangent = v_proj / (vx_i * dx / distance + vy_i * dy / distance)
        velocities[i] = (v_tangent * vx_i, v_tangent * vy_i)
    
    return velocities

if __name__ == "__main__":
    # Example usage
    # Assume positions list is already calculated from p1_position.py module
    initial_theta = 3.10 * 2 * np.pi  # Example initial theta
    positions = calculate_spiral_positions(initial_theta)  # Calculate node positions
    velocities = calculate_velocities(positions)  # Calculate node velocities

    # Print velocities
    for i, (vx, vy) in enumerate(velocities):
        if i % 50 == 0:
            print(f"Node {i + 1}: vx = {vx:.6f}, vy = {vy:.6f},"
                  f" |v| = {np.sqrt(vx**2 + vy**2):.6f} m/s")
