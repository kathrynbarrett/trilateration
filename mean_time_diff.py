'''
Calculate the mean and standard deviation of times in time_diff.txt.
'''

import numpy as np
import matplotlib.pyplot as plt

file = 'time_diff.txt'
data = np.genfromtxt(file)

time_diff_short = data[0]
time_diff_long = data[1]

mean_short = np.mean(time_diff_short)
mean_long = np.mean(time_diff_long)

std_short = np.std(time_diff_short)
std_long = np.std(time_diff_long)

print(f'Short method: mean time difference = {mean_short}, stdev = {std_short}')
print(f'Long method: mean time difference = {mean_long}, stdev = {std_long}')

labels = ['Short method', 'Longer method']
plt.figure()
plt.hist(data, label = labels)
plt.legend()
plt.show()