import numpy as np

def var_differential_longitudinal_motion(D, d, lamda, r0):
    '''
    From 'The ESO differential image motion monitor', Sarazin and Roddier, 1990
    
    Parameters
        D: aperture (m)
        d: separation of apertures (m)
        lamda: wavelength (m)
        r0: Fried parameter (m)        
    
    Returns
        variance: variance of differential longitudinal motion
    
    '''

    b = 0.179*D**(-1/3) - 0.0968*d**(-1/3)
    angle_var = 2*(lamda**2)*(r0**(-5/3))*b
    phase_var = angle_var/((2*lamda/(np.pi*D))**2)
    print(angle_var, phase_var)
    return angle_var, phase_var


def var_differential_transverse_motion(D, d, lamda, r0):
    '''
    From 'The ESO differential image motion monitor', Sarazin and Roddier, 1990
    
    Parameters
        D: aperture (m)
        d: separation of apertures (m)
        lamda: wavelength (m)
        r0: Fried parameter (m)        
    
    Returns
        variance: variance of differential transverse motion
    
    '''

    b = 0.179*D**(-1/3) - 0.145*d**(-1/3)
    angle_var = 2*(lamda**2)*(r0**(-5/3))*b
    phase_var = angle_var/((2*lamda/(np.pi*D))**2)
    print(angle_var, phase_var)

    return angle_var, phase_var

def rad_to_arcsec(angle):
    '''
    Converts an angle in radians to arcseconds
    '''
    return angle * (180/np.pi) * 3600

def arcsec_to_rad(angle):
    '''
    Converts an angle in arcseconds to radians
    '''
    return angle / 3600 * (np.pi/180)