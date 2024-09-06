from utils import *

time = 100  # Time in seconds
theta_s = find_theta_at_time(time)
print("Theta_s at time", time, "s is", theta_s)
print("The number of turns ar time", time, "s is", theta_s / (2 * np.pi))

x, y = calculate_position_at_time(time)
print(f"Position of the first node at time {time} seconds: x = {x:.6f}, y = {y:.6f}")

positions = calculate_spiral_positions(theta_s)

for i in range(len(positions)):
    if i % 50 == 0:
        print("Position", i, "is", positions[i])

