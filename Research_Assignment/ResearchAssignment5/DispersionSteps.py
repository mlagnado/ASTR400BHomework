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
for i in Time_list:
    print(i)
#print(Time_list)
#for i in Time_list:
#    print(Jacobi_Rad(i, 'M31', 'MW'))

def Time_at_snap(list, galaxy1, galaxy2):
    '''
    This function takes in a list of snap numbers and returns the Time and separation of two galaxies at those snap numbers
    input:

    output:
        Sep
        Time
    '''
    Galaxy1 = np.genfromtxt(f'Orbit_{galaxy1}.txt')
    Galaxy2 = np.genfromtxt(f'Orbit_{galaxy2}.txt')
    Separation = vector_dif(Galaxy1, Galaxy2)
    All_Time = Galaxy1[:,0] #Is the same for galaxy 1 or 2
    Time = []
    Sep = []
    for i in list:
        Time.append(All_Time[i])
        Sep.append(Separation[i,0])
#    print(Sep, Time)
    return Sep, Time


#Read in orbits for plotting
MW = np.genfromtxt('Orbit_MW.txt')
M31 = np.genfromtxt('Orbit_M31.txt')
Seps, Times = Time_at_snap(Time_list, 'MW', 'M31')
fig = plt.figure()
plt.plot(MW[:,0], vector_dif(MW,M31)[:,0], label='MW-M31', color='r')
plt.scatter(Times, Seps)
plt.xlabel('Time [Gyr]')
plt.ylabel('Seperation [kpc]')
plt.title('Seperation between MW and M31 over time')
#plt.semilogy()
plt.legend()
plt.savefig('MW-M31sep.png')
plt.show()


#On each of those snap numbers produce a plot that shows the dispersion vs the radius
for i in range(len(Time_list)):
    disp, disp2, rad = calculate_dispersion(Time_list[i], 'MW','M31')

    fig = plt.figure()
    plt.scatter(rad,disp, s=2, label='Dispersion of MW')
    plt.scatter(rad,disp2, s=2, label='Dispersion of MW2')
    plt.ylabel('Dispersion [km/s]')
    plt.xlabel('Radius [kpc]')
    plt.title(f'Dispersion vs Radius at snap number: {Time_list[i]}')
    plt.legend()
    plt.show()