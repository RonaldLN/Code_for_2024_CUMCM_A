import numpy as np
from scipy.optimize import fsolve

def calculate_spiral_positions(
        initial_theta, 
        pitch=0.55, 
        d1=3.41, 
        d=2.2, 
        num_nodes=223
    ):
    """
    Calculate the positions of nodes along an Archimedean spiral path.
    
    Args:
    - initial_theta_turns (float): Initial theta for the first node in turns.
    - pitch (float): Spiral pitch in meters per turn.
    - d1 (float): Distance between the first and second node in meters.
    - d (float): Distance between subsequent nodes in meters.
    - num_nodes (int): Total number of nodes to calculate.
    
    Returns:
    - positions (list of tuples): List of (x, y) positions of the nodes.
    """
    # Spiral equation: r = pitch / (2 * pi) * theta
    r_of_theta = lambda theta: pitch / (2 * np.pi) * theta

    # Function to find the next theta_i+1 given theta_i and r_i, considering the distance d
    def next_theta(theta_i, r_i, distance):
        """
        Solves for the next theta (angle) such that the distance between the current node
        and the next node on the spiral is equal to the given distance.
        """
        def equation(theta_next):
            r_next = r_of_theta(theta_next)  # Use the spiral equation to calculate r_next
            # Distance equation between two points in polar coordinates
            return np.sqrt(r_i**2 + r_next**2 - 2 * r_i * r_next * np.cos(theta_next - theta_i)) - distance
        
        # Initial guess for theta_next: small increment from theta_i
        theta_next_initial_guess = theta_i + distance / r_i
        theta_next_solution = fsolve(equation, theta_next_initial_guess)[0]
        return theta_next_solution

    # Initialize the first node's position
    theta_1 = initial_theta  # Initial theta for the first node (in radians)
    r_1 = r_of_theta(theta_1)  # Radial distance for the first node

    # Convert the initial position to Cartesian coordinates
    x_1 = r_1 * np.cos(theta_1)
    y_1 = r_1 * np.sin(theta_1)
    positions = [(x_1, y_1, theta_1)]

    # Calculate positions for subsequent nodes
    theta_current = theta_1  # Start with the initial theta
    r_current = r_1  # Start with the initial radius

    # Loop through the nodes
    for i in range(num_nodes):
        # Use d1 for the first distance and d for all subsequent distances
        distance = d1 if i == 0 else d
        
        # Calculate the next theta based on the current theta and radius
        theta_next = next_theta(theta_current, r_current, distance)
        r_next = r_of_theta(theta_next)  # Calculate the next radius using the lambda function

        # Convert to Cartesian coordinates
        x_next = r_next * np.cos(theta_next)
        y_next = r_next * np.sin(theta_next)
        
        # Append the new position
        positions.append((x_next, y_next, theta_next))
        
        # Update current theta and radius for the next iteration
        theta_current = theta_next
        r_current = r_next

    return positions


if __name__ == "__main__":
    # Example usage
    initial_theta_turns = 16  # Theta value for the first node in turns
    positions = calculate_spiral_positions(initial_theta_turns * 2 * np.pi)

    # Print positions or save them as needed
    for i, (x, y, _) in enumerate(positions):
        print(f"Node {i + 1}: x = {x:.6f}, y = {y:.6f}")
