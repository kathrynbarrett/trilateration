
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy
import numpy as np
import scipy.integrate
import scipy.interpolate
import scipy.signal

def fits_read(file,i):
    hdu = fits.open(file,ignore_missing_simple=True)
    data = hdu[i].data
    hdr = hdu[i].header
    return data,hdr

def convolution(array1,array2):
    '''
    Convolves array1 with array2

    Returns
        con: Convolution
    '''
    print('step 1')
    a_fft = np.fft.fft(array1)
    print('step 2')
    b_fft = np.fft.fft(array2)
    print('step 3')
    c_fft = a_fft*b_fft
    print('step 4')
    con = np.fft.fftshift(np.fft.ifft(c_fft).real)
    print('done')
    return con

def shift(array,shift):
    '''
    Shift pulse. Returns array and shifted array, with added zeros to make arrays the same length.
    
    Parameters
        shift = pixels to shift

    Returns
        array1: input array with zeros at end
        array2: shifted input array (zeros at start)
    
    '''
    print('Shifting pulse')

    pix_shift = shift
    array1 = np.pad(array,pix_shift)
    array2 = np.roll(array1,pix_shift)

    return array1, array2

def find_time_dif(a1,a2):
    '''
    Finds the time difference between 2 arrays by finding the cross correlation and extracting the time difference (distance of peak of cross-correlation from centre)

    Parameters
        a1,a2: Arrays to find cross-correlation

    Returns
        time_diff: Time difference between 2 arrays 
    '''
    plt.figure()
    plt.plot(a1)
    plt.plot(a2)
    plt.show()

    end = convolution(a1,a2)

    plt.figure()
    plt.plot(end, linestyle = 'none', marker = '.')
    plt.axvline(int(len(end)/2)-1)
    print(len(end)/2-1)
    plt.show()
    
    pix_diff = int((len(end))/2) -1 - np.argmax(end)

    print(int((len(end))/2))
    print(np.argmax(end))
    print(f'Time difference: {pix_diff} s')

    return end,pix_diff
    


def middle_index(arr):
    return int((len(arr)-1)/2)



############### TESTED METHOD ##########################################
sodium_profile = np.array([3,4,2,3,1])
pulse = np.array([10,0,1,5,0])

return_sig = scipy.signal.convolve(sodium_profile,np.flip(pulse))
print(f'Return signal: {return_sig}')





# pad and shift 1 array
shifted1, shifted2 = shift(return_sig,1)


return_sig_auto = scipy.signal.convolve(shifted1, np.flip(shifted1))
print(f'Auto-correlation of return signal: {return_sig_auto}')

plt.figure()
plt.plot(return_sig_auto)
plt.show()


# add noise to shifted and non-shifted arrays
mult = 2
noisy1 = np.random.poisson(shifted1*mult)/mult
noisy2 = np.random.poisson(shifted2*mult)/mult

# convolve SECOND ARRAY MUST BE FLIPPED TO GET THE TIME DELAY
noisy_conv = scipy.signal.convolve(noisy1, np.flip(noisy2))
offset = middle_index(noisy_conv) - np.argmax(noisy_conv)

# convolve this with the auto-correlation of the 'perfect' signal
final_conv1 = scipy.signal.convolve(noisy_conv,np.flip(return_sig_auto))
final_conv2 = scipy.signal.convolve(noisy_conv,return_sig_auto)
offset = middle_index(final_conv1) - np.argmax(final_conv1)

plt.figure()
plt.plot(noisy_conv, label = 'Convolution of noisy signals with offset')
plt.plot(return_sig_auto, label = 'Auto-correlation of perfect signal')
plt.legend()
plt.show()

plt.figure()
plt.plot(final_conv1, label = 'flipped')
plt.plot(final_conv2, label = 'not flipped')
plt.legend()
plt.show()

print(offset)
##############################################################################



#plt.figure()
#plt.plot(noisy_conv)
#plt.show()


exit()













number = 4
data_file = 'sodium_profile_data/na_profile_clusters_20190621.fits'
res = 1500

data, _ = fits_read(data_file,0)
profile = np.array([data[0,number,:],data[1,number,:]])

height = profile[0]
strength = profile[1]

# interpolate sodium profile for desried resolution
f = scipy.interpolate.interp1d(height, strength, bounds_error=False, fill_value=0)
height_new = np.arange(min(height), max(height), res)
strength_new = f(height_new)
print(strength_new)
profile_interp = np.array([height_new, strength_new])

# length of array - create pulse array to be same length
array_lengths = len(height_new)
print(f'Number of pixels in interpolated sodium layer: {array_lengths}')

plt.figure()
plt.plot(height_new,strength_new)
plt.show()

a1, a2 = shift(strength_new,10)

plt.figure()
plt.plot(np.convolve(a1,a2))
plt.show()


plt.figure()
plt.plot(a1)
plt.plot(a2)
plt.show()


find_time_dif(a1,np.flip(a2))