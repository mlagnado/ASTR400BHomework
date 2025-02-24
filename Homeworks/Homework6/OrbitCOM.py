

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def Orbit_COM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
          galaxy (string) name of the galaxy
	  start (int) first snapshot to be read in
	  end (int) last snapshot to be read in
	  n (int) intervals to return COM
    outputs: 
	  
    """
    
    # compose the filename for output
    fileout = f"Orbit_{galaxy}.txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec

    delta = 0.1
    volDec = 2.0
    if galaxy == 'M33':
        volDec = 4.0
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start,end,n)
    '''
    BONUS:
        check if the input is eligible
    '''
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        print(f"Currently computing snap number {snap_id}, counter: {i}")
        # compose the data filename (be careful about the folder)
        snap_num = '000' + str(snap_id) #Turns the snap number into a string
        snap_num = snap_num[-3:] #Makes sure the snap number is only 3 characters long
        filename = (galaxy) + "/" + "%s_"%(galaxy) + snap_num + '.txt' #Creates the filename
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_pos = COM.COM_P(delta,volDec)
        COM_vel = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2])
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i][0] = COM.time.value/1000 #units of Gyr
        orbit[i][1] = COM_pos[0].value #.value takes away units
        orbit[i][2] = COM_pos[1].value
        orbit[i][3] = COM_pos[2].value
        orbit[i][4] = COM_vel[0].value
        orbit[i][5] = COM_vel[1].value
        orbit[i][6] = COM_vel[2].value
        
        # print snap_id to see the progress
        print(f"Finishing calculation on snap number {snap_id}")
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
#Orbit_COM("MW",0,801,5) #Only do if want to redo calculation ~3min
#Orbit_COM("M31",0,801,5)
#Orbit_COM("M33",0,801,5)

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MW = np.genfromtxt('Orbit_MW.txt')
print('MW table', MW.shape)
print(MW)
M31 = np.genfromtxt('Orbit_M31.txt')
print('M31 table', M31.shape)
print(M31)
M33 = np.genfromtxt('Orbit_M33.txt')
print('M33 table', M33.shape)
print(M33)

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def vector_dif(vector1, vector2):
    '''
    This functin computes the difference of two vectors
    Inputs:
        vector1 (array) array of position and velocity coords output from Orbit_COM()
        vector2 (array of same shape as vector 1) vector to be subtracted from vector1 also output from Orbit_COM()
        PV (string) tells wether we want the position of velocity
    Outputs:
        magnitude (float) the magnitude of the difference of the two vectors
    '''
    magnitude = np.zeros([len(vector1),2])
    for i in range(len(vector1)):
        x = vector1[i][1]-vector2[i][1] #difference in x pos
        y = vector1[i][2]-vector2[i][2] #difference in y pos
        z = vector1[i][3]-vector2[i][3] #difference in z pos
        vx = vector1[i][4]-vector2[i][4] #difference in x vel
        vy = vector1[i][5]-vector2[i][5] #difference in y vel
        vz = vector1[i][6]-vector2[i][6] #difference in z vel
        magnitude[i][0] = np.sqrt(x**2+y**2+z**2)
        magnitude[i][1] = np.sqrt(vx**2+vy**2+vz**2)
    return magnitude


# Determine the magnitude of the relative position and velocities 

# of MW and M31
#print(vector_dif(MW,M31))
# of M33 and M31
#print(vector_dif(M33,M31))



# Plot the Orbit of the galaxies 
#################################
fig = plt.figure()
plt.plot(MW[:,0], vector_dif(MW,M31)[:,0])
plt.xlabel('Time [Gyr]')
plt.ylabel('Seperatino [kpc]')
plt.title('Seperation between MW and M31 over time')


# Plot the orbital velocities of the galaxies 
#################################

