import numpy as np
from ReadFile import Read
import astropy.units as u
from CenterOfMass import CenterOfMass

def ComponentMass(filename, particle_type, distance):
	'''
	This function sums the total mass of one specified particle type
	inputs: filename (string) where the data will be read from
		particle_type (int) the type of particle that you wish to sum
			Halo = 1, Disk = 2, Bulge = 3
	output: total_mass (astropy quantity) the total mass of a particular particle type in the data
	'''
	time, total, data = Read(filename) #reads the data
	COM = CenterOfMass(filename, 2)
	COM_pos = COM.COM_P(0.1,2.0) #Obtain COM coordinates
	data['x'] = data['x'] - COM_pos[0].value #Shifts coordinates to be centered on COM
	data['y'] = data['y'] - COM_pos[1].value
	data['z'] = data['z'] - COM_pos[2].value
	r = np.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
	index = np.where(r <= distance) #Only looking for mass within the distance from the galaxy
	data = data[index]
	mass = np.sum(data[data['type'] == particle_type]['m'])
	return mass*1e10

#Code used to get each value in table
#Changed file between MW, M31, and M33 and component from 1,2,and 3
MWHaloMass = ComponentMass('MW/MW_000.txt',1,30)
#print(MWHaloMass)
