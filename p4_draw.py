import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file to check its content
file_path = 'positions.csv'
data = pd.read_csv(file_path)

# Display the first few rows of the dataframe to understand its structure
data.head()

# The data has two columns, which seem to be coordinates for points.
# Let's rename the columns to 'X' and 'Y' for clarity and then plot the data.

# Rename the columns
data.columns = ['X', 'Y']

# Plot the points and connect them with lines
plt.figure(figsize=(10, 6))
plt.plot(data['X'], data['Y'], marker='o')
plt.title('Plot of Points with Lines Connecting Them')
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.grid(True)
plt.show()
