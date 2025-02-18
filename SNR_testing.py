import numpy as np
import matplotlib.pyplot as plt
import scipy

def find_time_dif(a1,a2):
    end = scipy.signal.correlate(a1,a2)
    time_diff = ((len(end)-1)/2 - np.argmax(end))
    print(f'Time difference: {time_diff} points')
    return time_diff

def shift(array,pix_shift):
    emp = np.zeros(pix_shift)
    array1 = np.append(array,emp)
    array2 = np.append(emp,array)
    return array1, array2

def add_shot_noise(array1,nphot):

        heights = array1[1:]
        total = np.sum(heights)

        # Convert height to photon intensity
        n_phot = nphot # total number of photons
        scaling_factor = n_phot / total
        intensity = heights*scaling_factor

        # check
        if int(sum(intensity)) !=n_phot:
            print('Integral is not equal to photon number')

        # add poisson noise
        noisy = np.random.poisson(intensity)

        plt.figure()
        plt.plot(noisy, label = 'shot noise')
        plt.plot(intensity, label = 'no noise')
        plt.legend()
        plt.show()

        return noisy

y = np.array([0,2,6,19,45,27,34,75,45,10,5,0])

arr1,arr2 = shift(y,1)
n = 3
noisy_arr1 = add_shot_noise(arr1,n)
noisy_arr2 = add_shot_noise(arr2,n)

find_time_dif(noisy_arr1,noisy_arr2)