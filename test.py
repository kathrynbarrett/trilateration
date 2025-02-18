import tools

h_eff = 2000
angle = tools.arcsec_to_rad(10)

aperture = 0.4
sep = angle * h_eff
print(f'Separation: {sep} m')
wavelength = 500e-9
fried_param = 0.1


angle, phase = tools.var_differential_longitudinal_motion(aperture, sep,wavelength, fried_param)

print(f"Angle stdev: {tools.rad_to_arcsec((angle)**0.5)} arcsec")