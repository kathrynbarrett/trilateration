import matplotlib.pyplot as plt
import numpy as np


def deltat_to_theta(deltat,baseline):
    rad = deltat * 3e8/(2*baseline)
    arcsec = rad * (180/np.pi) * 3600
    return arcsec

plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rcParams.update({
    'font.size': 16,  # Global font size
    'axes.titlesize': 18,  # Font size for titles
    'axes.labelsize': 18,  # Font size for axis labels
    'xtick.labelsize': 16,  # Font size for x-axis tick labels
    'ytick.labelsize': 16,  # Font size for y-axis tick labels
    'legend.fontsize': 16,  # Font size for the legend
})


#name = 'peak_uncertainty'
#name = 'testing'
name = 'quartic'

frequency = 100 #Hz

#file7 = f'{name}_1e-07.txt'
file8 = f'{name}_1e-08.txt'
#file9 = f'{name}_1e-09.txt'
file10 = f'{name}_1e-10.txt'
#file11 = f'{name}_1e-11.txt'
file12 = f'{name}_1e-12.txt'


#data7 = np.genfromtxt(file7)
data8 = np.genfromtxt(file8)
#data9 = np.genfromtxt(file9)
data10 = np.genfromtxt(file10)  
#data11 = np.genfromtxt(file11)  
data12 = np.genfromtxt(file12)


#data = {f'data{t_res}_all1e{exp}': [] for t_res in range(7, 13) for exp in range(2, 13)}
data = {f'data{t_res}_all1e{exp}': [] for t_res in [8,10,12] for exp in range(2, 13)}

data_arrays = {
    #7: data7,
    8: data8,
    #9: data9,
    10: data10,
    #11: data11,
    12: data12
}

#for t_res in range(7, 13): # loop over timing resolutions
for t_res in [8,10,12]: # loop over timing resolutions
    for exp in range(2, 11):    # loop over exponents
        key = f'data{t_res}_all1e{exp}' # key for the dictionary
        all = []
        for n in range(len(data_arrays[t_res])): # loop over the data to find data with the correct exponent
            if int(data_arrays[t_res][n, 1]) == 10**exp:
                all.append(data_arrays[t_res][n, 2])
        if all == []:
            continue
        else:  
            mean = np.mean(all)
            uncertainty = np.std(all)
            data[key].append((mean, uncertainty))


print(data)
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
colors =  plt.get_cmap('tab10').colors  # Use the tab10 colormap

plt.figure(figsize=(9, 6))
#for t_res in range(7, 13): # loop over timing resolutions
for t_res in [8,10,12]: # loop over timing resolutions
    exps = []
    means = []
    uncertainties = []
    for exp in range(2, 9):    # loop over exponents
        key = f'data{t_res}_all1e{exp}' # key for the dictionary
        if data[key] == []:
            print(key)
        else:
            mean = deltat_to_theta(data[key][0][0], baseline=10000)
            uncertainty = deltat_to_theta(data[key][0][1], baseline=10000)
            exps.append(10**exp * frequency)
            means.append(mean)
            uncertainties.append(uncertainty)
    print(f'exps: {exps}')
    print(f'means: {means}')
    print(f'uncertainties: {uncertainties}')
    # Plot the line without the label
    plt.errorbar(exps, means, yerr=uncertainties, fmt='-', capsize=6, linestyle='--', ms=5, color = colors[t_res-8],lw=2)
    # Plot the markers with the label
    plt.plot(exps, means, marker='o', label=f'1e-{t_res}', linestyle='None', ms=3, color = colors[t_res-8])  

plt.yscale('log')
plt.xscale('log')
plt.xlabel('Return photon rate / s$^{-1}$')
plt.ylabel('$\\theta_{\\mathrm{min}}$ / arcsec')
#plt.axhline(0.0003, color='black', linestyle='dotted')
plt.axvline(1e7, color='black', linestyle='dotted')
#plt.fill_between([0,1e7], 0, 10e7, color='gray', alpha=0.3)
plt.fill_between([10**3, 10**13], 0, 3e-5, color='black', alpha=0.9)
plt.xlim(5 * 10**3, 2 * 10**10)
plt.ylim(0, 0.2 *10e4)
plt.yticks([10**(-4), 10**(-2), 10**(0), 10**(2), 10**(4)])
plt.legend(title = 'Timing resolution / s', loc='upper right', bbox_to_anchor=(0.35, 0.4))
plt.grid()

plt.tick_params(axis='both', which='both', direction='out', length=6, width=0.5, colors='black', grid_color='gray', grid_alpha=0.5)
plt.savefig('peak_uncertainty.pdf', bbox_inches='tight')
plt.show()


exit()


x07 = data07[:,0] * frequency 
y07 = data07[:,1] * 3e8/2000 * (180/np.pi) * 3600 

x08 = data08[:,0] * frequency
y08 = data08[:,1] * 3e8/2000 * (180/np.pi) * 3600 

x09 = data09[:,0] * frequency
y09 = data09[:,1] * 3e8/2000 * (180/np.pi) * 3600 

x10 = data10[:,0] * frequency
y10 = data10[:,1] * 3e8/2000 * (180/np.pi) * 3600

#x11 = data11[:,0] * frequency
#y11 = data11[:,1] * 3e8/2000 * (180/np.pi) * 3600

#x12 = data12[:,0] * frequency
#y12 = data12[:,1] * 3e8/2000 * (180/np.pi) * 3600


plt.figure()
plt.plot(x07, y07, label='10$^{-7}$', marker='o', ms = 3)
plt.plot(x08, y08, label='10$^{-8}', marker='o', ms = 3)
plt.plot(x09, y09, label='1e-9', marker='o', ms = 3)
plt.plot(x10, y10, label='1e-10', marker='o', ms = 3)
#plt.plot(x11, y11, label='1e-11', marker='o', ms = 3)
#plt.plot(x12, y12, label='1e-12', marker='o', ms = 3)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Return photons / s')
plt.ylabel('$\\theta_{min}$ / \'\'')
plt.axvline(1e7, color='black', linestyle='dashed')
plt.axhline(0.0003, color='black', linestyle='dashed')
plt.legend(title = 'Timing resolution / s')
plt.grid()
plt.savefig('peak_uncertainty.pdf')
plt.show()