
# # Lab 5 ASTR 400B 
# 



# Import Modules 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const # import astropy constants
import astropy.units as u

print('## Part A ##')
# # Part A :  Mass to Light Ratios 
# 
# Wolf et al. 2010 
# 
# $M(<R_{half}) = \frac {4}{G}\sigma^2 R_e$
# 
# Where $R_{half}$ = 3D half mass radius 
# and $R_e$ is the 2D half mass radius of stars (observed)
# 
# Determine which of the following two systems are galaxies:
# 
# The system 47 Tuc is observed with:  $\sigma = 17.3$ km/s, $R_e = 0.5$ pc, $L_v \sim 10^5 L_\odot$ 
# 
# The system Willman I is observed with: $\sigma = 4.3$ km/s, $R_e = 25$ pc, $L_v = 10^3 L_\odot$



# Gravitational Constant in the desired units
# kpc^3/Gyr^2/Msun
Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)




def WolfMass(sigma, re):
    """ Function that defines the Wolf mass estimator from Wolf+ 2010
    PARAMETERS
    ----------
        sigma: astropy quantity
            1D line of sight velocity dispersion in km/s
        re: astropy quantity
            Effective radius, 2D radius enclosing half the
            stellar mass in kpc
    OUTPUTS
    -------
        mWolf: Returns the dynamical mass within the 
            half light radius in Msun
    """
    
    sigmaKpcGyr = sigma.to(u.kpc/u.Gyr) # velocity dispersion units
    
    mWolf = 4/Grav*sigmaKpcGyr**2*re # Wolf mass estimator
    
    return mWolf

##47 Tuc param
LumTuc = 1*10**5 * u.Lsun #luminosity
SigmaTuc = 17.1*u.km/u.s #1d loss velocity dispersion
reTuc = 0.5/1000*u.kpc #effective radius

#dynamical mass for 47 Tuc
massTuc = WolfMass(SigmaTuc, reTuc)
print(f"Dynamical mass for 47 Tuc: {massTuc:.2e}")

print(f"Mass to light ratio of 47 Tuc {np.round(massTuc/LumTuc,1)}")

LumWI = 1e3*u.Lsun #luminosity
sigmaWI = 4.3*u.km/u.s # los velocity dispersion
reWI = 25/1000*u.kpc # effective radius

#dynamical mass
massWI = WolfMass(sigmaWI, reWI)
print(f"Dynamical mass for Willman I: {massWI:.2e}")

print(f"Mass to light ratio of Willman I  {np.round(massWI/LumWI,1)}")
#This is missing a lot of mass (not close to 1) so we can tell it is a galaxy.




print('## Part B ##')
# # Part B :  Stellar to Halo Mass Relation
# 
# Following the work of [Moster et al. 2013 (MNRAS, 428, 3121)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.428.3121M/abstract)
# 
# 
# `Equation 2:`                  $ \frac{m}{M} = 2N \left [ \left ( \frac{M}{M_1} \right)^{-\beta} + \left (\frac{M}{M_1} \right)^{\gamma} \right]$ 
# 
# $m$ = stellar mass, $M$ = halo mass
# 
# `Equation 11:`        log $M_1(z) = M_{10} + M_{11} \frac{z}{z+1} $ 
# 
# `Equation 12:`        $N(z) = N_{10} + N_{11} \frac{z}{z+1} $
# 
# `Equation 13:`         $\beta(z) = \beta_{10} + \beta_{11} \frac{z}{z+1} $
# 
# `Equation 14:`         $\gamma(z) = \gamma_{10} + \gamma_{11} \frac{z}{z+1} $
print('## Question 1 ##')
# # Q1 
# 
# Modify the class below by adding a function called `StellarMass` that uses the `SHMratio` function and returns the stellar mass.



class AbundanceMatching:
    """ Class to define the abundance matching relations from 
    Moster et al. 2013, which relate the stellar mass of a galaxy
    to the expected dark matter halo mass, according to 
    Lambda Cold Dark Matter (LCDM) theory """
    
    
    def __init__(self, mhalo, z):
        """ Initialize the class
        
        PARAMETERS
        ----------
            mhalo: float
                Halo mass in Msun
            z: float
                redshift
        """
        
        #initializing the parameters:
        self.mhalo = mhalo # Halo Mass in Msun
        self.z = z  # Redshift
        
        
    def logM1(self):
        """eq. 11 of Moster 2013
        OUTPUT: 
            M1: float 
                characteristic mass in log(Msun)
        """
        M10      = 11.59
        M11      = 1.195 
        return M10 + M11*(self.z/(1+self.z))  
    
    
    def N(self):
        """eq. 12 of Moster 2013
        OUTPUT: 
            Normalization for eq. 2
        """
        N10      = 0.0351
        N11      = -0.0247
    
        return N10 + N11*(self.z/(1+self.z))
    
    
    def Beta(self):
        """eq. 13 of Moster 2013
        OUTPUT:  power of the low mass slope"""
        beta10      = 1.376
        beta11      = -0.826
    
        return beta10 + beta11*(self.z/(1+self.z))
    
    def Gamma(self):
        """eq. 14 of Moster 2013
        OUTPUT: power of the high mass slope """
        gamma10      = 0.608
        gamma11      = 0.329
    
        return gamma10 + gamma11*(self.z/(1+self.z))
    
    
    def SHMratio(self):
        """ 
        eq. 2 of Moster + 2013
        The ratio of the stellar mass to the halo mass
        
        OUTPUT: 
            SHMratio float
                Stellar mass to halo mass ratio
        """
        M1 = 10**self.logM1() # Converting characteristic mass 
        # to Msun from Log(Msun)
        
        A = (self.mhalo/M1)**(-self.Beta())  # Low mass end
        
        B = (self.mhalo/M1)**(self.Gamma())   # High mass end
        
        Norm = 2*self.N() # Normalization
    
        SHMratio = Norm*(A+B)**(-1)
    
        return SHMratio 
    
    def StellarMass(self):
        '''
        To compute the stellar mass using equation 2 of Moster + 2013
        (stellar/halo mass ratio)

        no inputs
        outputs:
          Starmass float stellar mass in Msun
        '''
        starmass = self.mhalo * self.SHMratio()
        return starmass

print('## Question 1 ##')
 # Q1: add a function to the class that takes the SHM ratio and returns 
# The stellar mass 


print('## Part C ##')
# # Part C : Plot the Moster Relation
# 
# Reproduce the below figure from Moster + 2013 
# Plot this for z=0, 0.5, 1, 2
# 
# ![mos](./MosterFig.png)



mh = np.logspace(10,15,1000) # Logarithmically spaced array
#halo mass array



# Define Instances of the Class for each redshift
MosterZ0 = AbundanceMatching(mh,0)
MosterZ05 = AbundanceMatching(mh,0.5)
MosterZ1 = AbundanceMatching(mh,1)
MosterZ2 = AbundanceMatching(mh,2)


fig,ax = plt.subplots(figsize=(10,8))


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Plot z = 0
plt.plot(np.log10(mh), np.log10(MosterZ0.StellarMass()),
         linewidth = 5, label='z=0')
plt.plot(np.log10(mh), np.log10(MosterZ05.StellarMass()),
         linewidth = 5, label='z=0.5',linestyle='--')
plt.plot(np.log10(mh), np.log10(MosterZ1.StellarMass()),
         linewidth = 5, label='z=1',linestyle=':')
plt.plot(np.log10(mh), np.log10(MosterZ2.StellarMass()),
         linewidth = 5, label='z=2',linestyle='-.')
# Continue plotting for the other redshifts here




# Axes labels 
plt.xlabel('log (M$_h$/M$_\odot$)',fontsize=22) 
plt.ylabel('log (m$_\star$/M$_\odot$)', fontsize=22)

# Legend
plt.legend(loc='lower right',fontsize='x-large')

# save the file 
plt.savefig('AbundanceMatching_Lab5.png')

print('## Part D ##')
# # Part D
# 
print('## Question 1 ##')
# # Q1
# 
# In studies that have modeled the Magellanic Clouds prior to 2010, the LMC is traditioanlly modeled with a halo (dark matter) mass of order $3 \times 10^{10}$M$_\odot$.  
# 
print('# A #')
# ## A) 
# According to $\Lambda$CDM theory, what should be the stellar mass of the LMC halo be at z=0?  
halo_LMC1 = 3*10**10 #traditional model
#create an abundance matching object
LMC1 = AbundanceMatching(halo_LMC1, 0)
#Find the stellar mass for this object
LMC1star = LMC1.StellarMass()
print(np.round(LMC1star/1e9,3))
print(np.round(LMC1star/3e9*100))
# 
print('# B #')
# ## B) 
# How does this stellar mass compare to the actual observed stellar mass of the LMC at the present day of ~$3 \times 10^9$ M$_\odot$ ? 
halo_LMC2 = 17e10
LMC2 = AbundanceMatching(halo_LMC2,0)
LMC2star = LMC2.StellarMass()
print(np.round(LMC2star/1e9,3))
# 
print('# C #')
# ## C) 
# What is the $\Lambda$CDM expected halo mass for the LMC (using Abundance Matching)?  





print('## Question 2 ##')
# # Q2
# 
# ## A) 
# What is the expected stellar mass of an L* galaxy at z=0? 
# 
M1halo_z0 = MosterZ0.logM1()
print(f"Log M1, z=0: {M1halo_z0}")

#Create a new instance of the class with halo mass = logM1 at z=0
M1Z0 = AbundanceMatching(10**M1halo_z0,0)
M1star_z0 = M1Z0.StellarMass()
print(f'Stellar mass of an L* galaxy at z=0: {np.round(M1star_z0/1e10,3)}')
# ## B)
# What is the expected stellar mass of an L* galaxy at z = 2? 
#Same thing for Z=2
M1halo_z2 = MosterZ2.logM1()
print(f"Log M1, z=2: {M1halo_z2}")
M1Z2 = AbundanceMatching(10**M1halo_z2,2)
M1star_z2 = M1Z2.StellarMass()
print(f'Stellar mass of an L* galaxy at z=2: {np.round(M1star_z2/1e10,3)}')


