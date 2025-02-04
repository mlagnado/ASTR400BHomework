MWDisk = ComponentMass(MW,2)
MWBulge = ComponentMass(MW,3)
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


MW = 'MW_000.txt'
M31 = 'M31_000.txt'
M33 = 'M33_000.txt'
import pandas as pd
MWHalo = ComponentMass(MW,1)
MWDisk = ComponentMass(MW,2)
MWBulge = ComponentMass(MW,3)
MWtotal = MWHalo + MWDisk + MWBulge
M31Halo = ComponentMass(M31,1)
M31Disk = ComponentMass(M31,2)
M31Bulge = ComponentMass(M31,3)
M31total = M31Halo + M31Disk + M31Bulge
M33Halo = ComponentMass(M33,1)
M33Disk = ComponentMass(M33,2)
M33Bulge = ComponentMass(M33,3)
M33total = M33Halo + M33Disk + M33Bulge
data = {'Galaxy Name': ['MW','M31','M33'], 'Halo Mass [$10^{12}M_{sun}$]':[MWHalo, M31Halo, M33Halo],
 'Disk Mass [$10^{12}M_{sun}$]': [MWDisk, M31Disk, M33Disk],
 'Bulge Mass [$10^{12}M_{sun}$]': [MWBulge, M31Bulge, M33Bulge],
 'Total [$10^{12}M_{sun}$]': [MWtotal, M31total, M33total],
} # 'f_bar':[]}
df = pd.DataFrame(data=data)
print(df)
