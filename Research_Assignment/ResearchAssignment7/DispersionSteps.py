'''
Velocity Dispersion: How does the velocity dispersion of each galaxy, in this showcase MW, evolve (bulge/disk) throughout the interaction (until they merge 6.5 Gyr)
I am going to look at how the velocity dispersion of MW change throughout its interaction with M31.
Velocity Dispersion is a quantity that is useful to display how well ordered the galaxy is or how disordered it becomes after an interaction.
'''

#Goal of this code
'''
This code is designed to create a plot for the interesting snap numbers.
The interesting snap numbers occur at local minima and maxima for separation between MW and M31.
'''

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

#Import previous functions
from ReadFile import Read
from OrbitCOM import vector_dif
from GalaxyMass import ComponentMass
from CenterOfMass import CenterOfMass
from VelocityDispersion import Find_Timesteps
from VelocityDispersion import calculate_dispersion

#Create a list of snap numbers
Time_list =  Find_Timesteps('MW','M31')
#Galaxy1 = np.genfromtxt(f'Orbit_MW.txt')
#for i in Time_list:
#    print(i)
#print(Time_list)
#for i in Time_list:
#    print(Jacobi_Rad(i, 'M31', 'MW'))

def Time_at_snap(snaps, galaxy1, galaxy2):
    '''
    This function takes in an array of snap numbers and returns the Time and separation of two galaxies at those snap numbers
    input:
        snaps: array of integers, Includes all the snap numbers
        galaxy1: string, the galaxy that the dispersion will be calculated for
        galaxy2: string, the other galaxy that impacts jacobi radius
    output:
        Sep: array of floats, Includes the separation between two galaxies at each snap number
        Time: array of floats, Includes the time in Gyr at each snap number
    '''
    Galaxy1 = np.genfromtxt(f'Orbit_{galaxy1}.txt') #Loading in the data
    Galaxy2 = np.genfromtxt(f'Orbit_{galaxy2}.txt')
    Separation = vector_dif(Galaxy1, Galaxy2) #Outputs the difference in position and velocity
    All_Time = Galaxy1[:,0] #Is the same for galaxy 1 or 2
    Time = np.zeros(len(snaps))
    Sep = np.zeros(len(snaps))
    for i in range(len(snaps)): #get the times and separations for these time steps
        Time[i] = All_Time[snaps[i]]
        Sep[i] = Separation[snaps[i],0]
#    print(Sep, Time)
    return Sep, Time


#Read in orbits for plotting
MW = np.genfromtxt('Orbit_MW.txt')
M31 = np.genfromtxt('Orbit_M31.txt')
Seps, Times = Time_at_snap(Time_list, 'MW', 'M31')
fig = plt.figure()
plt.plot(MW[:,0], vector_dif(MW,M31)[:,0], label='MW-M31', color='r')
plt.scatter(Times, Seps, label='Local Max/Min Separation')
plt.xlabel('Time [Gyr]')
plt.ylabel('Seperation [kpc]')
plt.title('Seperation between MW and M31 over time')
plt.legend()
plt.savefig('MLagnado_Separation.png')
plt.show()


#On each of those snap numbers produce a plot that shows the dispersion vs the radius
#fig = plt.figure()
#cmap = plt.get_cmap('cool') 
#color_val = np.linspace(0,1,len(Time_list))
#for i in range(len(Time_list)):
#    disp, disp2, rad = calculate_dispersion(Time_list[i], 'MW','M31')
#    color = cmap(color_val[i])
#    plt.plot(rad[1:],disp[1:], linewidth=2, alpha=0.5, label=f'Dispersion of MW at t={Times[i]}Gyr', color=color)
#    plt.plot(rad[1:],disp2[1:], linewidth=2, alpha=0.5, label=f'Dispersion of MW at t={Times[i]}Gyr', color=color)

#plt.ylabel('Dispersion [km/s]')
#plt.xlabel('Radius [kpc]')
#plt.title(f'Dispersion vs Radius')
#plt.legend()
#plt.savefig(f'MW_VelocityDispersion1.png')
#plt.show()

fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(10,6))
cmap = plt.get_cmap('cool') 
color_val = np.linspace(0,1,len(Time_list))
for i in range(len(Time_list)):
    disp, disp2, rad = calculate_dispersion(Time_list[i], 'MW','M31')
    color = cmap(color_val[i])
    ax[0].plot(rad[1:],disp[1:], linewidth=2, alpha=0.5, label=f'Dispersion of MW at t={np.round(Times[i],1)}Gyr', color=color)

ax[0].set_ylabel('Dispersion [km/s]')
ax[0].set_xlabel('Radius [kpc]')
ax[0].set_title(f'Dispersion vs Radius')
ax[0].legend()


#print quantified values
#disp0, disp02, rad = calculate_dispersion(Time_list[0], 'MW', 'M31')
#dispL, dispL2, rad = calculate_dispersion(Time_list[1], 'MW', 'M31')
#difference = dispL2-disp02
#print(difference[-4:])
#print(rad[-4:])


disp0, disp02, rad = calculate_dispersion(Time_list[0], 'MW', 'M31')
for i in range(len(Time_list)):
    if i == 0:
        continue
    disp, disp2, rad = calculate_dispersion(Time_list[i], 'MW','M31')
    difference = disp2-disp0
    print(difference[-4:])
    color = cmap(color_val[i])
    ax[1].plot(rad[1:],difference[1:], linewidth=2, alpha=0.5, label=f'Dispersion Change of MW at t={np.round(Times[i],1)}Gyr', color=color)
print(rad[-4:])
ax[1].set_ylabel('Dispersion [km/s]')
ax[1].set_xlabel('Radius [kpc]')
ax[1].set_title(f'Dispersion Change vs Radius')
ax[1].legend()
plt.savefig(f'MLagnado_Dispersion.png')
plt.show()