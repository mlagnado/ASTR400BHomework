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
Velocity Dispersion: How does the velocity dispersion of each galaxy, in this showcase MW, evolve (bulge/disk) throughout the interaction (until they merge 6.5 Gyr)
I am going to look at how the velocity dispersion of MW change throughout its interaction with M31.
Velocity Dispersion is a quantity that is useful to display how well ordered the galaxy is or how disordered it becomes after an interaction.
'''

#Outline of the process
'''
We want to find how the dispersion changes as a function of radius. To do this we have to select interesting points in time.
Interesting points we be when the two galaxies (MW and M31) are at a local max or min in separation (very close or very far).
Next we want to define up to what radius will we measure the dispersion.
From the center up to the jacobi radius (farthest distance still gravitationally bound).
Because we are looking at dispersion as a function of radius we have to choose a way to divide up the distance.
If we take spherical shells with a thickness ~2kpc we will have segments to calculate the dispersion inside of.
After calculating dispersion we can plot it to see how well it works.
'''

#Create a list of local minima and maxima no need to distinguish between min/max 
def Find_Timesteps(G1, G2):
    '''
    This function finds the local max and local min of seperation between two galaxy
    input:
        G1/G2: string, abreviation for galaxy names to be loaded in.
            Should already have a file from homework 6 loaded into this folder
    output:
        snap_nums: np.array, containes the snap numbers that represent the timestep which holds relevant value
    '''
    #Output from Homework 6 which will be used to calculate local min and local max seperation
    Galaxy1 = np.genfromtxt(f'Orbit_{G1}.txt')
    Galaxy2 = np.genfromtxt(f'Orbit_{G2}.txt')

    Separation = vector_dif(Galaxy1, Galaxy2) #This fxn spits out the distance and velocity separation
    Separation = Separation[:,0]
    times = Galaxy1[:,0]
#    print(Separation)
    index = np.where(times <= 6.50)[0]  #Get the indices where times <= 6.5 Gyr
    Separation = Separation[index] #Filter Separation by index
#    print(Separation)
    times = times[index] #Filter time by that index
#    print(Separation[1:-1], Separation[:-2], Separation[2:])
    maxima_indices = np.where((Separation[1:-1] > Separation[0:-2]) & (Separation[1:-1] > Separation[2:]))[0] + 1  # Local maxima add 1 because we shift the current by 1
#    print(maxima_indices)
    minima_indices = np.where((Separation[1:-1] < Separation[0:-2]) & (Separation[1:-1] < Separation[2:]))[0] + 1  # Local minima

    # Combine maxima and minima indices
    snap_indices = np.concatenate([maxima_indices, minima_indices])
    
    snap_nums = np.array([0] + list(snap_indices))
    snap_nums = np.sort(snap_nums)

    return snap_nums

#To calculate the Jacobi radius we need to know the separation between MW and M31 at each of our points of interest
#Jacobi radius is calculated by:
            #$r_j = R_0 (\frac{m}{2M}) ^ 1/3$           m = mass galaxy we are calculating for
                                                       #M = total mass
#This calculates the separation between them (R_0) in the jacobi radius eqn.
def Separation(snap_num, G1, G2):
    '''
    This function grabs the distance that two galaxies are separated by
    input:
        snap_num: int, the snap number where the separation will be retrieved
    output:
        Separation_time: float, the separation between the two galaxies at the desired snap number
    '''
    Galaxy1 = np.genfromtxt(f'Orbit_{G1}.txt')
    Galaxy2 = np.genfromtxt(f'Orbit_{G2}.txt')    
    Separation_all = vector_dif(Galaxy1, Galaxy2) #get distance and velocity separation

    Separation_time = Separation_all[snap_num][0] #Provides the separation of two galaxies at a given snap number
    return Separation_time


def Jacobi_Rad(snap, galaxy1, galaxy2):
    '''
    This function calculates the jacobi radius at a specified timestep for an input galaxy
    inputs:
        snap: int, The snap number that the jacobi radius will be calculated at
        galaxy1: string, the name of the galaxy that the jacobi radius will be calculated for ('MW, 'M31'...)
        galaxy2: string, the name of the galaxy that is bound to the first ('MW, 'M31'...)    
    output:
        J_rad: float, the Jacobi radius of galaxy1

    Eqn from King et. al. 1962
    '''
    snap_num = '000' + str(snap) #Turns the snap number into a string
    snap_num = snap_num[-3:] #Makes sure the snap number is only 3 characters long
    G1_filename = galaxy1 + '/' + galaxy1 + '_' + snap_num + '.txt' #opens the correct file
    G1_Halo = ComponentMass(G1_filename, 1) #calculates mass of component
    G1_Disk = ComponentMass(G1_filename, 2)
    G1_Bulge = ComponentMass(G1_filename, 3)
    G1_mass = G1_Bulge + G1_Disk + G1_Halo #adds mass of all components

    #same is done for the other galaxy
    G2_filename = galaxy2 + '/' + galaxy2 + '_' + snap_num + '.txt'
    G2_Halo = ComponentMass(G2_filename, 1)
    G2_Disk = ComponentMass(G2_filename, 2)
    G2_Bulge = ComponentMass(G2_filename, 3)
    G2_mass = G2_Bulge + G2_Disk + G2_Halo

    Tot_mass = G1_mass + G2_mass #gets the total mass for jacobi calc

    R_0 = Separation(snap, galaxy1, galaxy2) #Utilizes previous fxn for separation

    rj = R_0 * (G1_mass / (2 * Tot_mass))**(1/3) #As in lab 4 calculate jacobi radius
    return rj


#Now we calculate the dispersion a lot goes into this
#First we will get all the necessary coordinates and center them on MW
#Next we create a list of radii that contain the information for all the shells we will iterate over
#Finally we iterate over each shell and calculate velocity dispersion
#Dispersion is calculated by:
                        #sigma = sqrt(1/N * sum((v_i - V_mean)^2))
                        #V_mean is the mean velocity within the shell
def calculate_dispersion(snap,galaxy1, galaxy2,r=2.0/3):
    '''
    This function
    input:
        snap: int, the snap number that the calculation will occur at
        galaxy1: string, the galaxy that the dispersion will be calculated for
        galaxy2: string, the other galaxy that impacts jacobi radius
        r: float, softening length
    output:
        dispersion_list: np.array of floats, the dispersion within each of the shells
        dispersion_list2: same as above but calculated slightly differently
        radius_list: np.array of floats, the radii at the edge of each of the shells
    '''
    snap_num = '000' + str(snap) #Turns the snap number into a string
    snap_num = snap_num[-3:] #Makes sure the snap number is only 3 characters long
    filename = galaxy1 + '/' + galaxy1 + '_' + snap_num + '.txt'
    time, total, data = Read(filename)
    index = np.where(data['type'] != 1) #Find disk and bulge particles by looking at which ones are not halo
    COM = CenterOfMass(filename, 2)
    COM_pos = COM.COM_P(0.1,2.0) #Obtain COM coordinates
    COM_vel = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2]) #obtain COM velocity coordinates 
    
    m = data['m'][index] #only looking at dispersion of visible particles
    x_pos = data['x'][index]*u.kpc #get x position data for chosen particle type
    y_pos = data['y'][index]*u.kpc
    z_pos = data['z'][index]*u.kpc
    x_vel = data['vx'][index]*u.km/u.s
    y_vel = data['vy'][index]*u.km/u.s
    z_vel = data['vz'][index]*u.km/u.s


    x_new = x_pos - COM_pos[0] #Shifts coordinates to be centered on  galaxy 1 COM
    y_new = y_pos - COM_pos[1]
    z_new = z_pos - COM_pos[2]
    R = np.sqrt(x_new**2 + y_new**2 + z_new**2) #Now we have a list of all the particles in galaxy 1 and their distance to the center of that galaxy
    vx_new = x_vel - COM_vel[0] #Get velocity in terms of the COM
    vy_new = y_vel - COM_vel[1]
    vz_new = z_vel - COM_vel[2]
    vR = np.sqrt(vx_new**2 + vy_new**2 + vz_new**2)

    jacobi_radius = Jacobi_Rad(snap, galaxy1, galaxy2)#Calculates the jacobi radius at this timestep

    #Produce a list of radii from the center(0) until one step past the jacobi radius
    radii = np.arange(0.0,jacobi_radius+r,r) #at each of these radius we will calculate the dispersion inwards
                            #We add an r to make sure we don't get any errors later when iterating over radius.
#    print(jacobi_radius,radii)

#Calculates the dispersion in two ways.
    i = 1
    #We will go to the 'edge' of MW when the two galaxies are far apart
    upper_rad = 30
    if jacobi_radius < 30:
        upper_rad = jacobi_radius
    #First method: obtain magnitude of velocity and calculate dispersion of the magnitude of 3D velocity
    #Second method: obtains the magnitude of each velocity vector component and calculates the magnitude of the velocity dispersion from them
    dispersion_list = np.zeros(len(radii[radii < upper_rad])) #initializes lists for outputs
    dispersion2_list = np.zeros(len(radii[radii < upper_rad])) 
    radius_list = np.zeros(len(radii[radii < upper_rad])) 
    #We will first iterate over all radius starting from the edge of the first shell until the edge of the jacobi radius
    while radii[i] < upper_rad:
        current_radii = R.value #create a temporary list with R coordinates
        current_vel = vR.value #again temporary list of velocity coords
        current_vx = vx_new.value
        current_vy = vy_new.value
        current_vz = vz_new.value
#        print(radii[i-1], radii[i])
#        print(current_vel)
        index_low = np.where(current_radii >= radii[i-1]) #Finds an index for the lower edge of the shell
        current_radii = current_radii[index_low] #Gets rid of values from our list that are below this lower bound
        current_vel = current_vel[index_low]
        current_vx = current_vx[index_low]
        current_vy = current_vy[index_low]
        current_vz = current_vz[index_low]
#        print(radii[i-1], radii[i])
        index_high = np.where(current_radii < radii[i]) #Finds an index for the upper edge of the shell
        current_radii = current_radii[index_high] #Gets rid of values beyond the upper edge of shell
        current_vel = current_vel[index_high]
        current_vx = current_vx[index_high]
        current_vy = current_vy[index_high]
        current_vz = current_vz[index_high]
#        print(len(current_vel))
        #If our bounds didn't catch any particles this makes sure we don't have an error
        #It is mostly used far beyond the galaxy >30kpc for some rare cases
        if len(current_vel) < 3: 
            i += 1
            continue
        #Now find avg velocity in this shell
        mean_V = np.sum(current_vel)/len(current_vel) #avg vel
        mean_Vx = np.sum(current_vx)/len(current_vx) #x avg vel
        mean_Vy = np.sum(current_vy)/len(current_vy)
        mean_Vz = np.sum(current_vz)/len(current_vz)

        
        #And calculate dispersion in this shell from the mean
        sqr_dif = 0
        sqr_difx = 0
        sqr_dify = 0
        sqr_difz = 0

        #calculates the square difference for the values in our bounds
        for k in range(len(current_vel)):
            sqr_dif += (current_vel[k] - mean_V)**2
            sqr_difx += (current_vx[k] - mean_Vx)**2
            sqr_dify += (current_vy[k] - mean_Vy)**2
            sqr_difz += (current_vz[k] - mean_Vz)**2
#            print(sqr_dif)
        #Lastly we calculate the dispersion for each of our velocities
        dispersion = np.sqrt((1/len(current_vel)) * sqr_dif)
        dispersionx = np.sqrt((1/len(current_vx)) * sqr_difx)
        dispersiony = np.sqrt((1/len(current_vy)) * sqr_dify)
        dispersionz = np.sqrt((1/len(current_vz)) * sqr_difz)

        #Now for our second method we find the magnitude of this dispersion
        dispersion2 = np.sqrt(dispersionx**2 + dispersiony**2 + dispersionz**2)

        #Now add all of these values to the list
        dispersion_list[i] = dispersion
        dispersion2_list[i] = dispersion2
        radius_list[i] = radii[i]

        #continue iterating
        i += 1

    return dispersion_list, dispersion2_list, radius_list
'''
Testing it out
'''
#Create a list of snap numbers
#Time_list =  Find_Timesteps('MW','M31')
#Galaxy1 = np.genfromtxt(f'Orbit_MW.txt')
#for i in Time_list:
#    print(Galaxy1[i][0])
#print(Time_list)
#for i in Time_list:
#    print(Jacobi_Rad(i, 'M31', 'MW'))

#On each of those snap numbers produce a plot that shows the dispersion vs the radius
#for i in range(len(Time_list)):
#    disp, disp2, rad = calculate_dispersion(Time_list[i], 'MW','M31')

#    fig = plt.figure()
#    plt.scatter(rad,disp, s=2, label='Dispersion of MW')
#    plt.scatter(rad,disp2, s=3, label='Dispersion of MW2')
#    plt.ylabel('Dispersion [km/s]')
#    plt.xlabel('Radius [kpc]')
#    plt.title('Dispersion vs Radius')
#    plt.legend()
#    plt.show()

#Only want to save a plot of one of them
#disp, disp2, rad = calculate_dispersion(Time_list[0], 'MW','M31')
#fig = plt.figure()
#plt.scatter(rad,disp, s=3, label='Dispersion of MW')
#plt.scatter(rad,disp2, s=3, label='Dispersion of MW2')
#plt.ylabel('Dispersion [km/s]')
#plt.xlabel('Radius [kpc]')
#plt.title('Dispersion vs Radius')
#plt.legend()
#plt.savefig(f'MW_Dispersion_{Time_list[0]}')
#plt.show()