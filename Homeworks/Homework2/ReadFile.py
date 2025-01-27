import numpy as np
import astropy.units as u

def Read(file_name):
	file = open(file_name, 'r')
	line1 = file.readline()
	label, value = line1.split()
	time = float(value)*u.Myr
	line2 = file.readline()
	label2, value2 = line2.split()
	total = value2
	file.close()

	data = np.genfromtxt(file_name, dtype=None, names=True, skip_header=3)
	return time, total, data

#This is just a test of the function
#time, total, data = Read('MW_000.txt')
#print(data['y'][0])
