


# In Class Lab 6
# Surface Brightness Profiles




# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass


# # Lab 6: Sersic Profiles
# 
# In this lab we will use Homework 5 solutions to compute the mass profile of M31's bulge. 
# 
# We will turn the mass profile into a density profile and see if we can fit it reasonably well with a sersic profile. 

# # Part A
# 
# Determine the Surface Mass Density Profile for the Simulated M31 bulge
# 
# a) find the positions and masses of all the bulge particles



# Create a center of mass object for M31
# I.e. an instance of the CenterOfMass class 




# Use the center of mass object to 
# store the x, y, z, positions and mass of the bulge particles
# be sure to correct for the COM position of M31
 


# b) Change the bulge particle positions into cylindrical coordinates



# Determine the positions of the bulge particles in 
# cylindrical coordinates. 


# c) determine the projected surface mass density profile



def SurfaceDensity(r,m):
    """ Function that computes the surface mass density profile
    given an array of particle masses and radii 
     
    PARMETERS
    ---------
        r : array of `floats` - cyclindrical radius [kpc]
        m : array of `floats` - particle masses [1e10 Msun] 
    
    RETURNS
    -------
        r_annuli : array of `floats` -  radial bins for the 
            annuli that correspond to the surface mass density profile
    
        sigma: array of `floats` - surface mass density profile 
         [1e10 Msun/kpc^2] 
        
       
    """
    
    # Create an array of radii that captures the extent of the bulge
    # 95% of max range of bulge
    radii = np.arange(0.1, 0.95 * r.max(), 0.1)

    # create a mask to select particles within each radius
    # np.newaxis creates a virtual axis to make cyl_r_mag 2 dimensional
    # so that all radii can be compared simultaneously
    # a way of avoiding a loop - returns a boolean 
    enc_mask = r[:, np.newaxis] < radii

    # calculate mass of bulge particles within each annulus.  
    # relevant particles will be selected by enc_mask (i.e., *1)
    # outer particles will be ignored (i.e., *0)
    # axis =0 flattens to 1D
    m_enc = np.sum(m[:, np.newaxis] * enc_mask, axis=0)

    # use the difference between m_enc at adjacent radii 
    # to get mass in each annulus
    m_annuli = np.diff(m_enc) # one element less then m_enc
    
    
    # Surface mass density of stars in the annulus
    # mass in annulus / surface area of the annulus. 
    # This is in units of 1e10
    sigma = m_annuli / (np.pi * (radii[1:]**2 - radii[:-1]**2))
    # array starts at 0, but here starting at 1 and
    # subtracting radius that ends one index earlier.
    
    # Define the range of annuli
    # here we choose the geometric mean between adjacent radii
    r_annuli = np.sqrt(radii[1:] * radii[:-1]) 

    return r_annuli, sigma




# Define the surface mass density profile for the simulated bulge
# and the corresponding annuli


# # Part B 
# 
# 

# a) Create a function called `sersicE` that returns the Sersic Profile in terms of the effective radius $R_e$ (i.e. the half light radius).
# 
# $I(r) = I_e exp^{-7.67 ( (r/R_e)^{1/n} - 1)}$
# 
# Where 
# 
# $ L = 7.2 I_e \pi R_e^2$
# 
# and  $R_e$ is the half light radius.  We will assume a mass to light ratio for the stellar bulge of 1, so this is also the half mass radius.
# 
# The function should take as input the radius, $R_e$, $n$ and the total stellar mass of the system.
# 



def sersicE(r, re, n, mtot):
    """ Function that computes the Sersic Profile for an Elliptical 
    System, assuming M/L ~ 1. As such, this function is also the 
    mass surface density profile. 
    
    PARMETERS
    ---------
        r: `float`
            Distance from the center of the galaxy (kpc)
        re: `float`
            The Effective radius (2D radius that contains 
            half the light) (kpc)
        n:  `float`
            sersic index
        mtot: `float`
            the total stellar mass (Msun)

    RETURNS
    -------
        I: `array of floats`
            the surface brightness/mass density
            profile for an elliptical in Lsun/kpc^2

    """


    I=0
    
    return I
    
    


# b) Compute the M31 Bulge Mass Profile using Homework 5 MassProfile
# 



# Create a mass profile object for M31
# using solution to Homework 5




# Determine the Bulge mass profile
# use the annuli defined for the surface mass density profile


# c) Determine the effective radius of the bulge. Compute the total mass using Component Mass, from the GalaxyMass code, and find the radius that contains half this mass. 



# Determine the total mass of the bulge




# Find the effective radius of the bulge, 
# Re encloses half of the total bulge mass

# Half the total bulge mass




# Find the indices where the bulge mass is larger than b_half


# take first index where Bulge Mass > b_half
# check : should match b_half




# Define the Effective radius of the bulge


# d) Define the Sersic Profile for the M31 Bulge



# Sersic Index = 4


# # Part D
# 
# a) Plot the simulated surface mass density profile of M31's bulge as a proxy for its surface brightness profile (M/L = 1) 
# 
# b) Plot the M31 Bulge Sersic Profile
# 



fig, ax = plt.subplots(figsize=(9, 8))

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Surface Density Profile
# YOU ADD HERE


# Sersic fit to the surface brightness Sersic fit
# YOU ADD HERE



plt.xlabel('log r [ kpc]', fontsize=22)

# note the y axis units
plt.ylabel(r'log $\Sigma_{bulge}$ [$10^{10} M_\odot$ / kpc$^2$]', 
          fontsize=22)

plt.title('M31 Bulge', fontsize=22)

#set axis limits
plt.xlim(1,50)
plt.ylim(1e-5,0.1)

ax.legend(loc='best', fontsize=22)
fig.tight_layout()

#plt.savefig('Lab6.png')






