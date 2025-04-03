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
'''
'''



def Find_Timesteps(G1, G2):
    '''
    This function finds the local max and local min of seperation between two galaxy
    input:
        G1/G2: string, abreviation for galaxy names to be loaded in.
            Should already have a file from homework 6 loaded into this folder
    output:
        time_index: np.array, containes integers that represent the timestep which holds relevant value
    '''
    Galaxy1 = np.genfromtxt(f'Orbit_{G1}.txt')
    Galaxy2 = np.genfromtxt(f'Orbit_{G2}.txt')

    #Output from Homework 6 which will be used to calculate local min and local max seperation
    Separation = vector_dif(Galaxy1, Galaxy2)

    time_index = [0] #We will look at the first timestep for sure and then loop to add other important timesteps
    i = 1
    print(Galaxy1[i][0])
    while Galaxy1[i][0] <= 6.75: #6.75 Gyr is when the close encounters start to occur often this is "during the merger"
        if Separation[i][0] > Separation[i-1][0] and Separation[i][0] > Separation[i+1][0]:
            time_index.append(i)
        elif Separation[i][0] < Separation[i-1][0] and Separation[i][0] < Separation[i+1][0]:
            time_index.append(i)
        i += 1

    return np.array(time_index)


def Separation(time_index, G1, G2):
    '''
    
    '''
    Galaxy1 = np.genfromtxt(f'Orbit_{G1}.txt')
    Galaxy2 = np.genfromtxt(f'Orbit_{G2}.txt')    
    Separation_all = vector_dif(Galaxy1, Galaxy2)

    Separation_time = Separation_all[time_index][0] #Provides the separation of two galaxies at a given time
    return Separation_time


def Jacobi_Rad(time_index, galaxy1, galaxy2):
    '''
    This function calculates the jacobi radius at a specified timestep for an input galaxy
    inputs:
        time: int, The timestep that the jacobi radius will be calculated at
        galaxy1: string, the name of the galaxy that the jacobi radius will be calculated for ('MW, 'M31'...)
        galaxy1: string, the name of the galaxy that is bound to the first ('MW, 'M31'...)    
    output:
        J_rad: 
    '''
    snap_num = '000' + str(time_index) #Turns the snap number into a string
    snap_num = snap_num[-3:] #Makes sure the snap number is only 3 characters long
    G1_filename = galaxy1 + '/' + galaxy1 + '_' + snap_num + '.txt'
    G1_Halo = ComponentMass(G1_filename, 1)
    G1_Disk = ComponentMass(G1_filename, 2)
    G1_Bulge = ComponentMass(G1_filename, 3)
    G1_mass = G1_Bulge + G1_Disk + G1_Halo

    G2_filename = galaxy2 + '/' + galaxy2 + '_' + snap_num + '.txt'
    G2_Halo = ComponentMass(G2_filename, 1)
    G2_Disk = ComponentMass(G2_filename, 2)
    G2_Bulge = ComponentMass(G2_filename, 3)
    G2_mass = G2_Bulge + G2_Disk + G2_Halo

    Tot_mass = G1_mass + G2_mass

    R = Separation(time_index, galaxy1, galaxy2)

    rj = R * (G1_mass / (2 * Tot_mass))**(1/3) #As in lab 4
    return rj



def calculate_dispersion(time_index,galaxy1, galaxy2,r=2.0):
    '''
    This function
    input:

        r: float, softening length
    '''
    snap_num = '000' + str(time_index) #Turns the snap number into a string
    snap_num = snap_num[-3:] #Makes sure the snap number is only 3 characters long
    filename = galaxy1 + '/' + galaxy1 + '_' + snap_num + '.txt'
    time, total, data = Read(filename)
    index = np.where(data['type'] != 1) #Find disk and bulge particles by looking at which ones are not halo
    COM = CenterOfMass(filename, 2)
    COM_pos = COM.COM_P(0.1,2.0)
    COM_vel = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2])
    
    m = data['m'][index] #only looking at dispersion of visible particles
    x_pos = data['x'][index]*u.kpc #get x position data for chosen particle type
    y_pos = data['y'][index]*u.kpc
    z_pos = data['z'][index]*u.kpc
    x_vel = data['vx'][index]*u.km/u.s
    y_vel = data['vy'][index]*u.km/u.s
    z_vel = data['vz'][index]*u.km/u.s


    x_new = x_pos - COM_pos[0]
    y_new = y_pos - COM_pos[1]
    z_new = z_pos - COM_pos[2]
    R = np.sqrt(x_new**2 + y_new**2 + z_new**2) #Now we have a list of all the particles in galaxy 1 and their distance to the center of that galaxy
    vx_new = x_vel - COM_vel[0] #Get velocity in terms of the COM
    vy_new = y_vel - COM_vel[1]
    vz_new = z_vel - COM_vel[2]
    vR = np.sqrt(vx_new**2 + vy_new**2 + vz_new**2)

    jacobi_radius = Jacobi_Rad(time_index, galaxy1, galaxy2)
    radii = np.arange(0.0,jacobi_radius+2,r) #at each of these radius we will calculate the dispersion inwards
#    print(jacobi_radius,radii)

    dispersion_list = []
    dispersion2_list = []
    radius_list = []
    i = 1
    while radii[i] < jacobi_radius:
        current_radii = R.value
        current_vel = vR.value
        current_vx = vx_new.value
        current_vy = vy_new.value
        current_vz = vz_new.value
#        print(current_vel)
        index_low = np.where(current_radii >= radii[i-1])
        current_radii = current_radii[index_low]
        current_vel = current_vel[index_low]
        current_vx = current_vx[index_low]
        current_vy = current_vy[index_low]
        current_vz = current_vz[index_low]
#        print(radii[i-1], radii[i])
        index_high = np.where(current_radii < radii[i])
        current_radii = current_radii[index_high]
        current_vel = current_vel[index_high]
        current_vx = current_vx[index_high]
        current_vy = current_vy[index_high]
        current_vz = current_vz[index_high]
        if len(current_vel) == 0:
            i += 1
            continue
        #Now find avg velocity in this shell
        mean_V = np.sum(current_vel)/len(current_vel)
        mean_Vx = np.sum(current_vx)/len(current_vx)
        mean_Vy = np.sum(current_vy)/len(current_vy)
        mean_Vz = np.sum(current_vz)/len(current_vz)

        #And calculate dispersion in this shell from the mean
        sqr_dif = 0
        sqr_difx = 0
        sqr_dify = 0
        sqr_difz = 0

        for k in range(len(current_vel)):
            sqr_dif += (current_vel[k] - mean_V)**2
            sqr_difx += (current_vx[k] - mean_Vx)**2
            sqr_dify += (current_vy[k] - mean_Vy)**2
            sqr_difz += (current_vz[k] - mean_Vz)**2
#            print(sqr_dif)
        dispersion = np.sqrt((1/len(current_vel)) * sqr_dif)
        dispersionx = np.sqrt((1/len(current_vx)) * sqr_difx)
        dispersiony = np.sqrt((1/len(current_vy)) * sqr_dify)
        dispersionz = np.sqrt((1/len(current_vz)) * sqr_difz)
        dispersion2 = np.sqrt(dispersionx**2 + dispersiony**2 + dispersionz**2)
        dispersion_list.append(dispersion)
        dispersion2_list.append(dispersion2)
        radius_list.append(R[i].value)
        i += 1

    return dispersion_list, dispersion2_list, radius_list

Time_list =  Find_Timesteps('MW','M31')
Galaxy1 = np.genfromtxt(f'Orbit_MW.txt')
for i in Time_list:
    print(Galaxy1[i][0])
print(Time_list)
#for i in Time_list:
#    print(Jacobi_Rad(i, 'M31', 'MW'))
disp, disp2, rad = calculate_dispersion(Time_list[0], 'MW','M31')

fig = plt.figure()
plt.scatter(rad,disp, label='Dispersion of MW')
plt.scatter(rad,disp2, label='Dispersion of MW2')
plt.ylabel('Dispersion [km/s]')
#plt.semilogy()
plt.xlabel('Radius [kpc]')
plt.title('Dispersion vs Radius')
plt.legend()
plt.show()