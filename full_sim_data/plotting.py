import matplotlib.pyplot as plt
import numpy as np

plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rcParams.update({
    'font.size': 12,  # Global font size
    'axes.titlesize': 16,  # Font size for titles
    'axes.labelsize': 16,  # Font size for axis labels
    'xtick.labelsize': 14,  # Font size for x-axis tick labels
    'ytick.labelsize': 14,  # Font size for y-axis tick labels
    'legend.fontsize': 14,  # Font size for the legend
})


name = 'peak_uncertainty'
#name = 'testing'

frequency = 100 #Hz

file07 = f'{name}_1e-07.txt'
file08 = f'{name}_1e-08.txt'
file09 = f'{name}_1e-09.txt'
file10 = f'{name}_1e-10.txt'
file11 = f'{name}_1e-11.txt'
file12 = f'{name}_1e-12.txt'



data10 = np.genfromtxt(file10)  
data11 = np.genfromtxt(file11)  
data12 = np.genfromtxt(file12)
data07 = np.genfromtxt(file07)
data08 = np.genfromtxt(file08)
data09 = np.genfromtxt(file09)


x10 = data10[:,0] * frequency
y10 = data10[:,1] * 3e8/2000 * (180/np.pi) * 3600 * 1000

x11 = data11[:,0] * frequency
y11 = data11[:,1] * 3e8/2000 * (180/np.pi) * 3600 * 1000

x12 = data12[:,0] * frequency
y12 = data12[:,1] * 3e8/2000 * (180/np.pi) * 3600 * 1000

x07 = data07[:,0] * frequency 
y07 = data07[:,1] * 3e8/2000 * (180/np.pi) * 3600 * 1000

x08 = data08[:,0] * frequency
y08 = data08[:,1] * 3e8/2000 * (180/np.pi) * 3600 * 1000

x09 = data09[:,0] * frequency
y09 = data09[:,1] * 3e8/2000 * (180/np.pi) * 3600 * 1000

plt.figure()
plt.plot(x07, y07, label='1e-7', marker='o', ms = 3)
plt.plot(x08, y08, label='1e-8', marker='o', ms = 3)
plt.plot(x09, y09, label='1e-9', marker='o', ms = 3)
plt.plot(x10, y10, label='1e-10', marker='o', ms = 3)
plt.plot(x11, y11, label='1e-11', marker='o', ms = 3)
plt.plot(x12, y12, label='1e-12', marker='o', ms = 3)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Return photons / s')
plt.ylabel('$\sigma _{\\theta}$ / mas')
plt.axvline(1e7, color='black', linestyle='dashed')
plt.axhline(0.3, color='black', linestyle='dashed')
plt.legend(title = 'Timing resolution / s')
plt.grid()
plt.savefig('peak_uncertainty.pdf')
plt.show()