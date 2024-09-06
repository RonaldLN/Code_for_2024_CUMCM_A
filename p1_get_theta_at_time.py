import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

initial_theta = 16 * 2 * np.pi  # Initial theta for the first node in radians

def calculate_arc_length(theta, pitch=0.55):
    """
    Calculate the arc length of the Archimedean spiral from theta = 0 to theta.
    """
    integrand = lambda theta: np.sqrt((pitch / (2 * np.pi))**2 + (pitch / (2 * np.pi) * theta)**2)
    arc_length, _ = quad(integrand, theta, initial_theta)  # Integrate from theta to initial_theta (16*2*pi)
    return arc_length

def find_theta_for_distance(s, pitch=0.55):
    """
    Find the value of theta that corresponds to a specific arc length (distance) s.
    """
    # Define the function to find root for
    equation = lambda theta: calculate_arc_length(theta, pitch) - s
    # Use a numerical solver to find the root
    theta_solution = fsolve(equation, s)[0]  # Initial guess is s
    return theta_solution

def calculate_position_at_time(t, pitch=0.55):
    """
    Calculate the position of the first node after time t seconds given constant speed of 1 m/s.
    
    Args:
    - t (float): Time in seconds.
    - pitch (float): Spiral pitch in meters per turn.
    
    Returns:
    - (x, y) (tuple of floats): Cartesian coordinates of the first node at time t.
    """
    # Distance traveled after time t
    s = t  # Since speed is 1 m/s, distance = time in meters
    
    # Find theta corresponding to this distance
    theta_s = find_theta_for_distance(s, pitch)
    
    # Calculate the radial distance r at theta_s
    r_s = pitch / (2 * np.pi) * theta_s
    
    # Convert to Cartesian coordinates
    x_s = r_s * np.cos(theta_s)
    y_s = r_s * np.sin(theta_s)
    
    return x_s, y_s

def find_theta_at_time(t, pitch=0.55):
    """
    Calculate the angle theta_s of the first node after time t seconds, given a constant speed of 1 m/s.
    
    Args:
    - t (float): Time in seconds.
    - pitch (float): Spiral pitch in meters per turn.
    
    Returns:
    - theta_s (float): Angle in radians of the first node at time t.
    """
    # Distance traveled after time t
    s = t  # Since speed is 1 m/s, distance = time in meters
    
    # Find theta corresponding to this distance
    theta_s = find_theta_for_distance(s, pitch)
    
    return theta_s


if __name__ == "__main__":
    # Example usage
    time = 100  # Time in seconds
    # x, y = calculate_position_at_time(time)
    # print(f"Position of the first node at time {time} seconds: x = {x:.6f}, y = {y:.6f}")
    theta_s = find_theta_at_time(time)
    print(f"Theta of the first node at time {time} seconds: theta_s = {theta_s:.6f} radians")
