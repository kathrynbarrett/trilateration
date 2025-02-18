"""
Code to see effect of convolving the laser pulse with the sodium layer profile

- can we extract picosecond timings?

"""

from astropy.io import fits
import matplotlib.pyplot as plt
import scipy
import numpy as np
import scipy.integrate
import scipy.interpolate
import scipy.signal
from scipy.optimize import curve_fit

plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rcParams.update({
    'font.size': 12,  # Global font size
    'axes.titlesize': 16,  # Font size for titles
    'axes.labelsize': 16,  # Font size for axis labels
    'xtick.labelsize': 14,  # Font size for x-axis tick labels
    'ytick.labelsize': 14,  # Font size for y-axis tick labels
    'legend.fontsize': 14,  # Font size for the legend
})

def fits_read(file,i):
    hdu = fits.open(file,ignore_missing_simple=True)
    data = hdu[i].data
    hdr = hdu[i].header
    return data,hdr


class PulseReturn():
    def __init__(self,number,res,plotting=False):

        '''
        Set the initial parameters - resolution of the sodium profile in metres
        Then interpolates the profile to achieve this desired resolution.

        Parameters
            number: Farley sodium profile number 
            res: Desired resolution of sodium profile and pulse [metres]
        '''

        self.plotting = plotting
        data_file = 'sodium_profile_data/na_profile_clusters_20190621.fits'
        data, _ = fits_read(data_file,0)
        self.profile = np.array([data[0,number,:],data[1,number,:]])
        self.res = res


        height = self.profile[0]
        strength = self.profile[1]


        # interpolate sodium profile for desried resolution
        f = scipy.interpolate.interp1d(height, strength)
        height_new = np.arange(min(height), max(height), self.res)
        strength_new = f(height_new)
        self.profile_interp = np.array([height_new, strength_new])

        # length of array - create pulse array to be same length
        self.array_lengths = len(height_new)
        print(f'Number of pixels in interpolated sodium layer: {self.array_lengths}')
        


    
    def square_pulse(self,t):
        """
        Simulate square wave pulse with pulse length t, and total array length = sodium profile length

        Parameters
            t: Pulse length

        """
        wavelength = 3e8*t # width of the pulse [m]
        pulse_pixel_width = wavelength/self.res # width of the pulse [pixels]

        pulse = np.zeros(self.array_lengths) # make pulse array length equal to sodium profile array length
        pulse[int(self.array_lengths/2 - pulse_pixel_width/2):int(self.array_lengths/2 + pulse_pixel_width/2)] = 0.1 # make central values = 1

        d = np.arange(0,len(pulse)*self.res,self.res)
        self.pulse = np.array([d,pulse])


    def find_time_dif(self,a1,a2):
        '''
        Finds the time difference between 2 arrays by finding the cross correlation and extracting the time difference (distance of peak of cross-correlation from centre)

        Parameters
            a1,a2: Arrays to find cross-correlation

        Returns
            time_diff: Time difference between 2 arrays 

        '''


        conv = self.convolution(a1,np.flip(a2)) # flip second array

        time_diff = (int((len(conv)-1)/2) - np.argmax(conv)) * self.res/(3e8)

        print(f'time diff (pixels): {int((len(conv)-1)/2) - np.argmax(conv)}')
        print(f'Time difference: {time_diff} s')

        return conv,time_diff
    
    
    def shift(self,array,time_difference):
        '''
        Shift pulse. Returns array and shifted array, with added zeros to make arrays the same length.
        
        Parameters
            time_difference = time to shift [s]

        Returns
            array1: input array with zeros at end
            array2: shifted input array (zeros at start)
        
        '''
        print('Shifting pulse')

        pix_shift = int(3e8*time_difference/self.res) 

        array1 = np.pad(array,pix_shift)
        array2 = np.roll(array1,pix_shift)

        #plt.figure()
        #plt.plot(array1)
        #plt.plot(array2)
        #plt.show()

        return array1, array2


    def convolution(self,array1,array2):

        '''
        Convolves array1 with array2

        Returns
            con: Convolution
        '''
        #conv = np.convolve(array1,array2,mode='same')

        a_fft = np.fft.fft(array1)
        b_fft = np.fft.fft(array2)
        c_fft = a_fft*b_fft
        con = np.fft.fftshift(np.fft.ifft(c_fft).real)

        return con

    
    
    def scale(self,rel_intensity,nphot):

        total = np.sum(rel_intensity)

        print('Converting return signal to intensity')

        # Convert height to photon intensity
        scaling_factor = nphot / total
        intensity = rel_intensity*scaling_factor

        # check
        if int(sum(intensity)) !=nphot:
            print('Integral is not equal to photon number')

        return intensity
    
    
    def add_shot_noise(self,intensity):
        print('Adding shot noise to return signal')
        # add poisson noise
        noisy = np.random.poisson(intensity)

        for vals in range(len(noisy)):
            if noisy[vals] <0:
                noisy[vals] = 0
        
        
        # for plotting, number 5
        distance1 = np.array(range(len(noisy)))
        time1 = distance1 / 3e8 * 1e6 
        #t = open(f'Data/return_sig_noisy2.txt', "w")
        #for i in range(len(c)):
        #    t.write(f'{time[i]} {noisy[i]}\n')
        #t.close()

        #if self.plotting:
        #    plt.figure()
        #    plt.plot(time1, noisy,label = 'shot noise', marker = '.', linestyle='none')
        #    plt.yticks([])
        #    plt.xlabel('Time')
        #    plt.ylabel('Flux')
        #    plt.title('Noisy return signal')
        #    plt.legend()
        #    plt.show()
        

        return noisy
    
    
    

##########################################################################################
num = 4 #Sodium profile number to use     
min_time_res = 1e-8 # resolution / s
resolution = min_time_res * 3e8 # resolution [metres/pixel] SHOULD BE 0.0003 for 1ps/pix
print(f'RESOLUTION = {resolution} m / pixel')
time = 1e-7 # pulse width [s]
tdiff = 0 #input time difference [s]
#photon_num = [10**4]
#photon_num = [10**3,10**4,10**5,10**6,10**7,10**8,10**9,10**10,10**11,10**12,10**13] # number of photons in pulse
photon_num = [10**6,10**7,10**8,10**9,10**10] # number of photons in pulse
##########################################################################################

'''
Steps:

1. Initialise class PulseReturn, with chosen sodium profile number and desired resolution in m/pix
2. Create laser pulse with given resolution

4. Convolve the pulse with the sodium profile
5. Scale so that the area under the curve is equal to the number of photons when running at picosecond frame rates 
    - decreases array size by 1
6. Shift one array to simulate time delay 
    - increases array size by number of pixels that correspond to the time shift
7. Find the autocorrelation of one of these arrays 
    - This assumes we know the exact pulse shape (ok) and the exact sodium profile (?)
    - This gives a 
8. Add shot noise to both arrays
9. Find the cross correlation
10. Find the difference between the cross-correlation and the auto-correlation using another auto-correlation

'''

#method = 'covariance'
#method = 'stdev'
method = 'quartic'

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))


def quartic(x, a, c, e, w):
    return a*(x-w)**4 + c*(x-w)**2 + e


if method == 'stdev':
    # Running 1000 times for a each photon number and timing resolution

    # initialise class
    a = PulseReturn(num,resolution,plotting=False)
    a.square_pulse(time)
    return_pulse = a.convolution(a.profile_interp[1],np.flip(a.pulse[1]))   




    # loop over photon numbers
    for its in photon_num:
        f = open(f'testing_{min_time_res}.txt', "a")

        # create array to store mean values
        mean_values = []

        intensity = a.scale(return_pulse,its)
        arr1,arr2 = a.shift(intensity,tdiff)



        # loop over 1000 times
        for i in range(1000):

            noisy_arr1 = a.add_shot_noise(arr1)
            noisy_arr2 = a.add_shot_noise(arr2)


            cross = a.convolution(noisy_arr1,np.flip(noisy_arr2))
            cross_times = np.array(range(len(cross)))

            cross = cross[int(2*len(cross)/5):int(3*len(cross)/5)]
            cross_times = cross_times[int(2*len(cross_times)/5):int(3*len(cross_times)/5)]

            cross = cross-min(cross)
            cross=cross/max(cross)  # rescale for gaussian fitting  

            popt, pcov = curve_fit(gaussian, cross_times, cross, p0=[1.0, np.mean(cross_times), np.std(cross)])

            # Extract the parameters and uncertainties
            amp, mean, sigma = popt
            mean_uncertainty = np.sqrt(pcov[1, 1])  # Uncertainty in the mean

            true_mean = (len(cross)/2) -1

            print(f"Fitted parameters:")
            print(f"Amplitude: {amp:.3f}")
            print(f"Mean (peak): {mean* resolution / 3e8} ± {mean_uncertainty* resolution / 3e8}")
            print(f"Standard deviation: {sigma* resolution / 3e8 }")

            plt.figure()
            plt.plot(cross_times*resolution/3e8,cross/max(cross), label = 'Data')  
            plt.plot(cross_times*resolution/3e8, gaussian(cross_times, *popt), color='red', lw=2, label="Fitted Gaussian")
            plt.errorbar(mean* resolution / 3e8, amp, xerr= mean_uncertainty*resolution/3e8, capsize=2, marker = 'x', color = 'black', label = 'Peak')
            plt.xlabel('Time / s')
            plt.ylabel('Convolution intensity')
            plt.grid()
            plt.legend()
            plt.savefig('fitted_gaussian.pdf', bbox_inches='tight')
            plt.show()

            mean_values.append((mean-true_mean)* resolution / 3e8)

        print(f'LENGTH of MEAN VALUES = {len(mean_values)}')


        print(np.mean(mean_values))
        plt.figure()
        plt.hist(mean_values,bins=50)
        
        plt.show()


        std_means = np.std(mean_values)
        f.write(f'{its} {std_means}\n')
        f.close()

elif method == 'covariance':
    # Running for different photon numbers for a given timing resolution

    #Loop through different sodium profiles
    for num in range(16):
        print(f'USING SODIUM PROFILE NUMBER {num}')

        # 1
        print('Interpolating sodium profile')
        a = PulseReturn(num,resolution,plotting=False)
        # 2
        print("Making square pulse")
        a.square_pulse(time) 
        # 4
        return_pulse = a.convolution(a.profile_interp[1],np.flip(a.pulse[1]))




        for its in photon_num:

            f = open(f'peak_uncertainty_{min_time_res}.txt', "a")
            # 5
            intensity = a.scale(return_pulse,its)


            # 6
            arr1,arr2 = a.shift(intensity,tdiff)

            # 8

            for i in range(len(arr1)):
                if arr1[i] < 0:
                    arr1[i] = 0

            for j in range(len(arr2)):
                if arr2[j] < 0:
                    arr2[j] = 0


            noisy_arr1 = a.add_shot_noise(arr1)
            noisy_arr2 = a.add_shot_noise(arr2)

            # 9
            print('Finding cross correlation between the 2 return signals')
            cross = a.convolution(noisy_arr1,np.flip(noisy_arr2))
            cross_times = np.array(range(len(cross))) #* resolution / 3e8

            #cross = cross[int(2*len(cross)/5):int(3*len(cross)/5)]
            #cross_times = cross_times[int(2*len(cross_times)/5):int(3*len(cross_times)/5)]

            #cross = cross-min(cross)

            cross=cross/max(cross)  # rescale for gaussian fitting  
            popt, pcov = curve_fit(gaussian, cross_times, cross, p0=[1.0, np.mean(cross_times), np.std(cross)])

            # Extract the parameters and uncertainties
            amp, mean, sigma = popt
            mean_uncertainty = np.sqrt(pcov[1, 1])  # Uncertainty in the mean


            print(f"Fitted parameters:")
            print(f"Amplitude: {amp:.3f}")
            print(f"Mean (peak): {mean* resolution / 3e8} ± {mean_uncertainty* resolution / 3e8}")
            print(f"Standard deviation: {sigma* resolution / 3e8 }")


            #plt.figure(figsize=(8, 5))
            #plt.plot(cross_times,cross, alpha=0.6, color='gray', label="Data")
            #plt.plot(cross_times, gaussian(cross_times, *popt), color='red', lw=2, label="Fitted Gaussian")
            #plt.xlabel("Value")
            #plt.ylabel("Density")
            #plt.legend()
            #plt.title("Gaussian Fit to Data")
            #plt.title(f'{its}')
            #plt.show()

        
            f.write(f'{num} {its} {mean_uncertainty* resolution / 3e8}\n')
            f.close()

elif method == 'quartic':
    # Running 1000 times for a each photon number and timing resolution

    for num in range(16):
        print(f'USING SODIUM PROFILE NUMBER {num}')

        # 1
        print('Interpolating sodium profile')
        a = PulseReturn(num,resolution,plotting=False)
        # 2
        print("Making square pulse")
        a.square_pulse(time) 
        # 4
        return_pulse = a.convolution(a.profile_interp[1],np.flip(a.pulse[1]))


        # loop over photon numbers
        for its in photon_num:
            f = open(f'quartic_{min_time_res}.txt', "a")

            # create array to store mean values
            mean_values = []

            intensity = a.scale(return_pulse,its)
            arr1,arr2 = a.shift(intensity,tdiff)

            for i in range(len(arr1)):
                if arr1[i] < 0:
                    arr1[i] = 0

            for j in range(len(arr2)):
                if arr2[j] < 0:
                    arr2[j] = 0

            noisy_arr1 = a.add_shot_noise(arr1)
            noisy_arr2 = a.add_shot_noise(arr2)


            cross = a.convolution(noisy_arr1,np.flip(noisy_arr2))
            cross_times = np.array(range(len(cross)))

            max_cross_times = max(cross_times)
            cross_times = cross_times/max_cross_times # normalise

            cross=cross/max(cross)
            #cross = cross[int(2*len(cross)/5):int(3*len(cross)/5)]
            #cross_times = cross_times[int(2*len(cross_times)/5):int(3*len(cross_times)/5)]


            # Fit the quartic polynomial
            initial_guess = [1, -0.5, 1, 0.5]  # Rough initial guess
            popt, pcov = curve_fit(quartic, cross_times, cross, p0=initial_guess)

            # Extract the fitted parameters
            coeffs = popt

            centre = coeffs[3]  # Centre of the peak
            centre_uncertainty = np.sqrt(np.abs(pcov[3, 3]))

            # Generate fitted values
            y_fit = quartic(cross_times, *popt)

            # Plot the original data and the fitted curve
            plt.scatter(cross_times, cross, label='Data', color='red', s=10)
            plt.plot(cross_times, y_fit, label='Quartic Fit', color='blue')
            plt.axvline(centre, color='green', linestyle='--', label=f'Peak: {centre:.3f} ± {centre_uncertainty:.3f}')
            plt.legend()
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.title('Quartic Fit with Peak Estimation')
            plt.show()

            peakx_rescaled = (centre) * max_cross_times * resolution /3e8
            uncertainty_on_peak_x_rescaled = centre_uncertainty * max_cross_times * resolution /3e8

            # Print results
            print(f'Peak at x = {peakx_rescaled:.3f} ± {uncertainty_on_peak_x_rescaled}')


        
            f.write(f'{num} {its} {uncertainty_on_peak_x_rescaled}\n')
            f.close()


'''


    #pixels_from_centre = (len(cross)-1)/2 - np.argmax(cross)
    #print(f'Pixels from centre: {pixels_from_centre}')
    #dis_diff = pixels_from_centre * resolution
    #print(f'Distance difference: {dis_diff}')
    #time_diff_short =dis_diff/(3e8)
    #print(f'Time difference: {time_diff_short} s')
    #f.write(f'{time_diff_short} ')


    #plt.figure()
    #plt.plot(cross,label='Cross-covariance')
    #plt.plot(auto,label='Auto-convariance')
    #plt.legend()
    #plt.xlabel('Time')
    #plt.xticks([])
    #plt.yticks([])
    #plt.title('Auto and cross')
    #plt.savefig('Auto-cross.pdf', bbox_inches='tight')
    #plt.show()

    # 10
    #print('Final step: Finding cross-correlation between auto and cross')
    #end, time_diff = a.find_time_dif(cross,auto)

    #end_times = np.array(range(len(end))) #* resolution / 3e8 

    #plt.figure()
    #plt.plot(end_times,end/max(end))
    #plt.show()
    #print(end_times[101])


    #end=end/max(end)  # rescale for gaussian fitting  
    #popt, pcov = curve_fit(gaussian, end_times, end, p0=[1.0, np.mean(end), np.std(end)])

    # Extract the parameters and uncertainties
    #amp, mean, sigma = popt
    #mean_uncertainty = np.sqrt(pcov[1, 1])  # Uncertainty in the mean


    #print(f"Fitted parameters:")
    #print(f"Amplitude: {amp:.3f}")
    #print(f"Mean (peak): {mean* resolution / 3e8} ± {mean_uncertainty* resolution / 3e8}")
    #print(f"Standard deviation: {sigma* resolution / 3e8 }")

    #plt.figure(figsize=(8, 5))
    #plt.plot(end_times,end, alpha=0.6, color='gray', label="Data")
    #plt.plot(end_times, gaussian(end_times, *popt), color='red', lw=2, label="Fitted Gaussian")
    #plt.xlabel("Value")
    #plt.ylabel("Density")
    #plt.legend()
    #plt.title("Gaussian Fit to Data")
    #plt.show()




    #plt.figure()
    #plt.title('Final')
    #plt.plot(a.convolution(cross,auto), linestyle = 'none', marker = '.')
    #plt.plot(a.convolution(cross,np.flip(auto)), linestyle = 'none', marker ='.')
    #plt.axvline(int((len(end)-1)/2))
    #plt.show()

'''
# for plotting, number 8
#distance = np.array(range(len(end)))
#time = distance/3e8 *1e6

#t = open(f'Data/cross_cor8.txt', "w")
#for i in range(len(time)):
#    t.write(f'{time[i]} {end[i]}\n')
#t.close()


#plt.figure()
#plt.plot(time,end)
#plt.title('Final autocorrelation')
#plt.show()
