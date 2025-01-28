from ReadFile import Read
import numpy as np
import astropy.units as u
import astropy.constants as const

def ParticleInfo(filename, particle_type, particle_number):
	time, total, data = Read(filename)
	if particle_type.lower() == 'dark matter':
		type = 1
	elif particle_type.lower() == 'disk stars':
		type = 2
	elif particle_type.lower() == 'bulge stars':
		type = 3
	index = np.where(data['type'] == type)
	print(index)
	xnew = data['x'][index]
	ynew = data['y'][index]
	znew = data['z'][index]
	pos_x = xnew[particle_number]
	pos_y = ynew[particle_number]
	pos_z = znew[particle_number]
	vxnew = data['vx'][index]
	vynew = data['vy'][index]
	vznew = data['vz'][index]
	vel_x = vxnew[particle_number]
	vel_y = vynew[particle_number]
	vel_z = vznew[particle_number]
	dist_mag = np.sqrt(pos_x**2 + pos_y**2 + pos_z**2)*u.kiloparsec
	vel_mag = np.sqrt(vel_x**2 + vel_y**2 + vel_z**2)*u.kilometer/u.second
	dist_mag = np.around(dist_mag,3)
	vel_mag = np.around(vel_mag,3)
	mass = data['m'][particle_number]*const.M_sun

	return dist_mag, vel_mag, mass

test_dist, test_vel, test_mass = ParticleInfo('MW_000.txt', 'Disk Stars', 100)
test_dist = np.around(test_dist.to(u.lyr),3)
print(test_dist, test_vel, test_mass)
