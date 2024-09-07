import numpy as np

def theta_to_spiral_position(theta, pitch=1.7):
    """
    Convert a given theta value to the (x, y) position on the spiral.

    Args:
        theta (float): The angle in radians.
        pitch (float): The pitch of the spiral.

    Returns:
        (float, float): The (x, y) position on the spiral.
    """
    r = pitch / (2 * np.pi) * theta
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y