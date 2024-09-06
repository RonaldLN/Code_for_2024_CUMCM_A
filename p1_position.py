import numpy as np
from scipy.optimize import fsolve

# Constants
pitch = 0.55  # Spiral pitch
d = 2.2  # Distance between nodes (in meters)

# Function to find the next theta_i+1 given theta_i and r_i
def next_theta(theta_i, r_i):
    def equation(theta_next):
        r_next = pitch * theta_next
        return np.sqrt(r_i**2 + r_next**2 - 2 * r_i * r_next * np.cos((theta_next - theta_i) * 2 * np.pi)) - d
    
    # Initial guess for theta_next could be theta_i + small increment
    theta_next_initial_guess = theta_i + d / (r_i * 2 * np.pi)
    theta_next_solution = fsolve(equation, theta_next_initial_guess)[0]
    return theta_next_solution

# Initialize the first node position
theta_1 = 16  # initial theta in radians
r_1 = pitch * theta_1

# Store the positions
positions = [(r_1 * np.cos(theta_1 * 2 * np.pi), r_1 * np.sin(theta_1 * 2 * np.pi))]

# Calculate subsequent nodes
num_nodes = 223  # Example for 223 nodes
theta_next = 0
for i in range(1, num_nodes):
    theta_i = theta_1 if i == 1 else theta_next
    r_i = pitch * theta_i
    
    # Calculate next theta
    theta_next = next_theta(theta_i, r_i)
    r_next = pitch * theta_next
    
    # Convert to Cartesian coordinates
    x_next = r_next * np.cos(theta_next * 2 * np.pi)
    y_next = r_next * np.sin(theta_next * 2 * np.pi)
    
    # Append the new position
    positions.append((x_next, y_next))

# Print positions or save them as needed
print(positions)
