import numpy as np
from ReadFile import Read
import astropy.units as u

def ComponentMass(filename, particle_type):
	'''
	This function sums the total mass of one specified particle type
	inputs: filename (string) where the data will be read from
		particle_type (int) the type of particle that you wish to sum
			Halo = 1, Disk = 2, Bulge = 3
	output: total_mass (astropy quantity) the total mass of a particular particle type in the data
	'''
	time, total, data = Read(filename) #reads the data
	mass = np.sum(data[data['type'] == particle_type]['m'])
	return np.round(mass*1e10/1e12, 3)

#Code used to get each value in table
#Changed file between MW, M31, and M33 and component from 1,2,and 3
MWHaloMass = ComponentMass('MW/MW_000.txt',1)
print(MWHaloMass)
