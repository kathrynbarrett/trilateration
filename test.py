import tools

h_eff = 4590
angle = tools.arcsec_to_rad(10)

aperture = 0.4
sep = angle * h_eff
print(f'Separation: {sep} m')
wavelength = 1550e-9
fried_param = 0.151


angle_long, phase_long = tools.var_differential_longitudinal_motion(aperture, sep,wavelength, fried_param)

angle_trans, phase_trans = tools.var_differential_transverse_motion(aperture, sep,wavelength, fried_param)

print(f"Longitudinal angle stdev: {tools.rad_to_arcsec((angle_long)**0.5)} arcsec")
print(f'Transverse angle stdev: {tools.rad_to_arcsec((angle_trans)**0.5)} arcsec')