import numpy as np
import matplotlib.pyplot as plt

# Example 1D data
y = np.array([1, 16, 81, 256, 625, 1296, 2401, 500, 0, -100, 50, 2000])  # Quartic data (example: y = x^4)
x = np.arange(len(y))  # Generate x values as indices

# Fit a quartic polynomial (degree 4)
coefficients = np.polyfit(x, y, 4)
polynomial = np.poly1d(coefficients)

# Generate fitted values
x_fit = np.linspace(x.min(), x.max(), 100)
y_fit = polynomial(x_fit)

# Plot the original data and the fitted curve
plt.scatter(x, y, label='Original Data', color='red')
plt.plot(x_fit, y_fit, label='Quartic Fit', color='blue')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Quartic Fit')
plt.show()

# Print the quartic equation
print(f'Fitted quartic equation: {polynomial}')
