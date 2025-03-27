# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 

# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
#from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, output_filename):
        '''
        This initializes the class and sets up various constant that will be used in the class calculations
        Calculates positions and velocities of M31 and M33
        Inputs:
            output_filename: string, the filename where the output data will be stored
        '''
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### store the output file name
        self.fileout = output_filename
        
        ### get the current pos/vel of M33 
        # create an instance of the  CenterOfMass class for M33 
        COM_M33 = CenterOfMass('M33_000.txt', 2)
        # store the position VECTOR of the M33 COM (.value to get rid of units)
        COM_M33_pos = COM_M33.COM_P(0.1,4.0) #in HW 6 we defiend volDec to be 4.0 for M33
        # store the velocity VECTOR of the M33 COM (.value to get rid of units)
        COM_M33_vel = COM_M33.COM_V(COM_M33_pos[0], COM_M33_pos[1], COM_M33_pos[2]).value
        COM_M33_pos = COM_M33_pos.value #needed units in COM_vel calculation so we remove it here
        ### get the current pos/vel of M31 
        # create an instance of the  CenterOfMass class for M31 
        COM_M31 = CenterOfMass('M31_000.txt',2)
        # store the position VECTOR of the M31 COM (.value to get rid of units)
        COM_M31_pos = COM_M31.COM_P(0.1,2.0)
        # store the velocity VECTOR of the M31 COM (.value to get rid of units)
        COM_M31_vel = COM_M31.COM_V(COM_M31_pos[0], COM_M31_pos[1], COM_M31_pos[2]).value
        COM_M31_pos = COM_M31_pos.value #needed units in COM_vel calculation so we remove it here
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0=np.array(COM_M33_pos-COM_M31_pos)
        self.v0=np.array(COM_M33_vel-COM_M31_vel)
        ### get the mass of each component in M31 
        ### disk
        # self.rdisk = scale length (no units)
        self.rdisk = 5 #*u.kpc
        # self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass('M31_000.txt',2)*1e12
        ### bulge
        # self.rbulge = set scale length (no units)
        self.rbulge = 1 #*u.kpc
        # self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass('M31_000.txt',3)*1e12
        # Halo
        # self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62 #*u.kpc
        # self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = 1.921*1e12
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        '''
        This function computes the Hernquist acceleration
        input:
            M: float, Mass of the galaxy component
            r_a: float, hernquist scale length
            r: np.array, position vector
        output:
            Hern: np.array, acceleration vector induced by a hernquist profile
        '''
        ### Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        
        ### Store the Acceleration
        Hern =  -self.G*M/(rmag*(r_a+rmag)**2) * r #because r is a np.array it does the calculation on each component
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r): # it is easiest if you take as an input a position VECTOR  r 
        '''
        This function calculates the acceleration using the Miyamoto Nagai profile (1975)
            This profile mimics the exponential disk far away from the disk
        input:
            M: float, Mass of galaxy component
            r_d: float, scale length for the disk of M31
            r: np.array, position vector
        output:
            Miya: np.array, acceleration vector
        '''
        
        ### Acceleration
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        # multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        R = np.sqrt(r[0]**2 + r[1]**2) #calculate R
        z_d = r_d/5.0 #z_d is defined as scale length/5
        B = r_d + np.sqrt(r[2]**2 +z_d**2) #calculate B
        Miya = -self.G*M/((R**2 + B**2)**1.5) * r #once again r is a np.array so calculation is on each component
        Miya *= np.array([1,1,B/np.sqrt(r[2]**2 + z_d**2)]) #z component is slightly more complicated so we have additional term to multiply

        return Miya
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r):
        '''
        This function sums up the acceleration due to every component
        input:
            r: np.array, position vector
        output:
            sum_accel: np.array, acceleration components summed
        '''
        ### Call the previous functions for the halo, bulge and disk
        # these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        Bulge_a = self.HernquistAccel(self.Mbulge, self.rbulge, self.r0) #calculate acceleration of the bulge
        Halo_a = self.HernquistAccel(self.Mhalo, self.rhalo, self.r0) #calculate acceleration of the halo
        Disk_a = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, self.r0) #calculate acceleration of the disk
        sum_accel = Bulge_a+Halo_a+Disk_a #sum the acceleration
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return sum_accel
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        '''
        This Function takes short timesteps to calculate r and v
        input:
            dt: float, size of each timestep
            r: np.array of floats, position vector
            v: np.array of floats, velocity vector
        output:
            rnew: np.array of floats, new position vector
            vnew: np.array of floats, new velocity vector
        '''
        if dt <= 0: #check that we are trying to move forward and not backwards in time
            print('please use a positive timestep value \n If you do want a negative timestep,')
            print('change the code in the LeapFrog function in the M33AnalyticOrbit class')
            rnew, vnew = r, v
            print('values are not being updated')
            return rnew, vnew
        # predict the position at the next half timestep
        rhalf = r + v*(dt/2)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = r + (1/2)*(v+vnew)*dt #update rnew
        
        return rnew, vnew #return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        '''
        This function performs an integration to calculate how M33 will move
        input:
            t0: float, initial time
            dt: float, size of each timestep
            tmax: float, end time
        '''
        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros([int(tmax/dt)+2, 7])
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        self.r = self.r0 #set current r
        self.v = self.v0 #set current/starting v
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t <= tmax): # as long as t has not exceeded the maximal time 
            
            # advance the time by one timestep, dt
            t += dt
            # store the new time in the first column of the ith row
            orbit[i][0] = t
            
            # advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            self.r, self.v = self.LeapFrog(dt, self.r, self.v) #calculate new r and v
         
    
            # store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i, 1:4] = self.r #store r
            
            # store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i, 4:8] = self.v #store current v
            
            # update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        
        # write the data to a file
        np.savetxt(self.fileout, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

#Calculat orbit
#Orbit = M33AnalyticOrbit('M33Orbit.txt')
#Orbit.OrbitIntegration(0,0.1,10)
M33_M31 = np.genfromtxt('M33Orbit.txt')
#create magnitude vector
def mag(vector):
    magnitude = np.zeros([len(vector),2])
    for i in range(len(vector)):
        magnitude[i][0] = np.sqrt(vector[i][1]**2+vector[i][2]**2+vector[i][3]**2)
        magnitude[i][1] = np.sqrt(vector[i][4]**2+vector[i][5]**2+vector[i][6]**2)
    return magnitude

OrbitMag = mag(M33_M31)
#import HW6 data
from OrbitCOM import vector_dif
M31 = np.genfromtxt('Orbit_M31.txt')
M33 = np.genfromtxt('Orbit_M33.txt')



#Plotting the orbit
fig, axes = plt.subplots(2, 1, figsize=(5, 8), sharex=True)
ax = axes[0]
ax.plot(M33_M31[:,0],OrbitMag[:,0], label='M33-M31 Analytic Orbit')
ax.plot(M31[:,0], vector_dif(M31,M33)[:,0], label='M33-M31 COM Orbit', color='r')
ax.semilogy()
ax.set(ylabel='Seperation [kpc]')
ax.set_title('Seperation between M33 and M31')
ax.legend()
ax=axes[1]
ax.plot(M33_M31[:,0],OrbitMag[:,1], label='M33-M31 Analytic Orbit')
ax.plot(M31[:,0], vector_dif(M31,M33)[:,1], label='M33-M31 COM Orbit', color='r')
ax.semilogy()
ax.set(xlabel='Time [Gyr]', ylabel='Relative Velocity [km/s]')
ax.set_title('Relative Velocity between M33 and M31')
ax.legend()
plt.savefig('M33-M31_Orbit')
plt.show()

print('## Question 2 ##')
print('Initially the plots are fairly similar however after the first inflection point the two graphs diverge.')
print('The main difference is that the Analytic Orbit only gives positive acceleration.')

print('## Question 3 ##')
print('In this problem it seems that we are accounting for the effects of M31 on M33 however not the effects of M33 on M31')
print('M33 should also be pulling M31 closer which means that they should be accelerating away from each other slower')
print('as opposed to M33 just slingshotting past a stationary M31.')

print('## Question 4 ##')
print('By considering the effects of MW, M33 would feel an acceleration in another direction as well. This means instead of falling right into M31,')
print('M33 would not get as close to it and would perhaps circle around it. I imagine that the Analytic Orbit plot for seperation might be wider')
print('and less sharp at the minima.')
