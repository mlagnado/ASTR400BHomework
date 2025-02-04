import numpy as np
from ReadFile import Read
import astropy.units as u

def ComponentMass(filename, particle_type):
	'''
	This function sums the total mass of one specified particle type
	inputs: filename (string) where the data will be read from
		particle_type (int) the type of particle that you wish to sum
			Halo = 1, Disk = 2, Bulge = 3
	output: total_mass (astropy value) the total mass of a particular particle type in the data
	'''
	time, total, data = Read(filename) #reads the data
	index = np.where(data['type'] == particle_type) #finds the indexes of one particle type
	mass_data = data['m'][index] #creates a new data list of mass from a single particle type
	total_mass = np.sum(mass_data)*u.M_sun #sums all the mass of a desired particle type
	total_mass = total_mass/(10**2) #mass is given in units of e10 MSol divide by e2 to get e12 units
	total_mass = np.round(total_mass, 3) #round the mass value to 3 decimal places
	return total_mass
