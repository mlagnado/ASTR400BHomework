#This functions calculates the 3D distance and 3D velocity of a specified particle
#it takes a file name as a string, particle type as a string, and specific particle as an int number as inputs
#the function returns: the 3D distance in kpc, 3D velocity in km/s, and mass in solar masses of the particle
from ReadFile import Read
import numpy as np
import astropy.units as u
import astropy.constants as const

def ParticleInfo(filename, particle_type, particle_number):
	time, total, data = Read(filename) #reads the file into the time, total number of particles, and a data array
	# the if statements try to convert the particle type into an integer so that they work with the data in the array
	if particle_type.lower() == 'dark matter':
		type = 1 #if the input is dark matter it will be type 1
	elif particle_type.lower() == 'disk stars':
		type = 2 #if the input is a disk star it will be type 2
	elif particle_type.lower() == 'bulge stars':
		type = 3 #if the input is a bulge star it will be type 3
	index = np.where(data['type'] == type) #creates an index to only look at particles of a specific type
	xnew = data['x'][index] #creates an array of the x positions of the chosen particle type
	ynew = data['y'][index] #creates an array of the y positions of the chosen particle type
	znew = data['z'][index] #creates an array of the z positions of the chosen particle type
	pos_x = xnew[particle_number] #takes the x position of the specified particle
	pos_y = ynew[particle_number] #takes the y position of the specified particle
	pos_z = znew[particle_number] #takes the z position of the specified particle
	vxnew = data['vx'][index] #creates an array of the x velocities of the chosen particle type
	vynew = data['vy'][index] #creates an array of the y velocities of the chosen particle type
	vznew = data['vz'][index] #creates an array of the z velocities of the chosen particle type
	vel_x = vxnew[particle_number] #takes the x velocity of the specified particle
	vel_y = vynew[particle_number] #takes the y velocity of the specified particle
	vel_z = vznew[particle_number] #takes the z velocity of the specified particle
	dist_mag = np.sqrt(pos_x**2 + pos_y**2 + pos_z**2)*u.kiloparsec #calculates the 3D distance of the particle and converts it to Kpc
	vel_mag = np.sqrt(vel_x**2 + vel_y**2 + vel_z**2)*u.kilometer/u.second #calculates the 3D velocity of the particle and converts the units to km/s
	dist_mag = np.around(dist_mag,3) #rounds the distance of the particle to 3 places
	vel_mag = np.around(vel_mag,3) #rounds the velocity of the particle to 3 places
	mnew = data['m'][index] #finds the mass of the particle and converts it to units of solar mass
	mass = mnew[particle_number]*(10**10)*u.kg
	mass = mass.to(u.M_sun)
	return dist_mag, vel_mag, mass

#testing the function works
test_dist, test_vel, test_mass = ParticleInfo('MW_000.txt', 'Disk Stars', 100)
test_dist = np.around(test_dist.to(u.lyr),3)
print('Answers to question 5')
print('distance:',test_dist)
print('velocity:', test_vel)
print('mass:', test_mass)
