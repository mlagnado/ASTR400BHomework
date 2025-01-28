from ReadFile import Read
import numpy as np
import astropy.units as u
import astropy.constants as const

def ParticleInfo(filename, particle_type, particle_number):
	time, total, data = Read(filename)
	if particle_type.lower() == 'dark matter'
		type = 1
	elif particle_type.lower() == 'disk stars'
		type = 2
	elif particle_type.lower() == 'bulge stars'
		type = 3
	index = np.where(data['#type'] = type)
	#come back later to do particle type part
	pos_x = data['x'][particle_number+index]
	pos_y = data['y'][particle_number+index]
	pos_z = data['z'][particle_number+index]
	vel_x = data['vx'][particle_number+index]
	vel_y = data['vy'][particle_number+index]
	vel_z = data['vz'][particle_number+index]
	dist_mag = np.sqrt(pos_x**2 + pos_y**2 + pos_z**2)*u.kiloparsec
	vel_mag = np.sqrt(vel_x**2 + vel_y**2 + vel_z**2)*u.kilometer/u.second
	dist_mag = np.around(dist_mag,3)
	vel_mag = np.around(vel_mag,3)
	mass = data['m'][particle_number]*const.M_sun

	return dist_mag, vel_mag, mass
