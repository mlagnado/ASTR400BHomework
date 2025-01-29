#This function reads a file to extract the data 
#The input is only a file name as a string
#The output are the  time, total number of particles and all the rest as a data frame
#The data frame includes particle type, mass, x,y,z distances, and vx,vy,vz velocities
import numpy as np
import astropy.units as u

def Read(file_name):
	file = open(file_name, 'r') #opens the file
	line1 = file.readline() #reads only the first line of the file
	label, value = line1.split() #splits the first line into two segments the title and value
	time = float(value)*u.Myr #takes the value from line one and converts it to a float with units (Myr)
	line2 = file.readline() #reads only the second line
	label2, value2 = line2.split() #splits the second line into two segments. A title and value
	total = value2 # assigns the value of the second line to the variable total
	file.close()

	data = np.genfromtxt(file_name, dtype=None, names=True, skip_header=3) #reads the rest of the file  into a data array
	return time, total, data

#This is just a test of the function
#time, total, data = Read('MW_000.txt')
#print(data['y'][0])
