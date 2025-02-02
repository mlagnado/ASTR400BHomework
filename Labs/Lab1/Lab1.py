
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" by midnight thursday

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants
# ### a)
# 
# Create a function called VLSR to compute the local standard of rest (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
RoReid = 8.34*u.kpc
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
RoAbuter = 8.178*u.kpc
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
RoSparke = 7.9*u.kpc
#
def VLSR(Ro, mu=6.379, vo=12.24*u.km/u.s):
    '''
    This function computes the velocity at the local standard of rest.
    VLSR = 4.74*mu*Ro-vsun
    Inputs Ro units of kpc distance form the sun to the galactic center
        mu is proper motion of Sag A* mas/yr default is from Reid & Brunthaler 2004
        vsun in units of km/s is the peculiar motion of the sun in the v direction from Schonrich 2010
    outputs: VLSR in units of km/s the local standard of rest
    '''
    vlsr = 4.74*mu*(Ro/u.kpc)*u.km/u.s - vo #computes the local standard of rest
    return vlsr

#1.
VLSR_Reid = VLSR(RoReid) #calculates the speed at the Reid estimate of the sun distance
#2.
VLSR_Abuter = VLSR(RoAbuter) #calculates the speed at the Abuter estimate of the sun distance
#3.
VLSR_Sparke = VLSR(RoSparke) #calculate the speed at the Sparke estimate of the sun distance
#just a test to see if it works
print('VLSR velocity Reid:',np.round(VLSR_Reid))
print('VLSR velocity Abuster:',np.round(VLSR_Abuter))
print('VLSR velocity Sparke:',np.round(VLSR_Sparke))




# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr
#V=2*pi*r/T
#T=2*pi*r/v
def TorbSun(Ro, Vc):
    '''
    This function calculates the orbital period of a particle at the distance of the sun
    inputs: radius (float) in kpc, speed (float) in km/s
    output orbital period (float) in Gyr in the "v" direction
    '''
    orbSun = (2*np.pi*Ro) / (Vc.to(u.kpc/u.Gyr)) #does the actual calculation converting the velocity to the correct units as well
    return orbSun
VsunPec = 12.24*u.km/u.s
Vsun = VLSR_Abuter + VsunPec

T_Abuter = TorbSun(RoAbuter,Vsun) #calculates the orbital period of a particle moving in a circular orbit at the suns position and velocity
print('Orbital period of the sun:', np.round(T_Abuter,3))



# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)
#1 rotation = 1 orbital period
#13.8Gyr/(1orbitalperiod/rotatio) = #rotations

rotations = np.round((13.8*u.Gyr)/T_Abuter,3) #calculates the number of rotations the sun could make about the GC  in 13.8 Gyr
print('Number of rotations:',rotations)



# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# point mass ~ 0
# star mass ~ 85
# dark matter ~ 11
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 
Grav = const.G.to(u.kpc**3/u.Gyr**2/u.M_sun)
#print(Grav)

#Density profile rho = VLSR**2/(4*pi*G*R**2)
#Mass (r) = Integrate rhodV
#       integrate rho 4*pi*r^2 dr
#       integrate VLSR^2 / (4*pi*G*R^2) * 4*pi*r^2 dr
#       VLSR^2/G *r

def massIso(r, VLSR):
    '''
    This function computes the dark matter mass enclosed by some radius (r) assuming an isothermal sphere model
    M(r) = VLSR^2/G*r
    inputs: r (astropy quantity) distance from galactic center (Kpc)
            VLSR (astropy quantity) velocity at the local standard of rest (Km/s)
    output: M (astropy quantity) mass enclosed within r (Msun)
    '''
    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) # converts units from km/s to kpc/gyr
    M = VLSRkpcGyr**2 /Grav *r #Isothermal sphere mass profile

    return M

#Compute the mass enclosed within the sun continue using Abuter values
MIsoSolar = massIso(RoAbuter, VLSR_Abuter)
print(f'Mass of the Isothermal model within Abuter Radius: {np.round(MIsoSolar,3):.2e}')

#mass enclosed by 260kpc
Miso260 = massIso(260*u.kpc, VLSR_Abuter)
print(f'Mass of the Isothermal model within 260 kpc: {np.round(Miso260,3):.2e}')
# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

#potential for a hernquist sphere
#phi = -G*M/ (r+a)

#escape speed becomes:
#Vesc^2 * 2+g+M/(r+a)

#rearrange for m
# M = vesc^2/2/G*(r+a)

def MassHernVesc(vesc,r,a=30*u.kpc):
    '''
    This functino determines the total dark mattter mass needed given an escape speed assuming a hernquist profile
    M = vesc^2/2/G*(r+a)
    inputs: vesc (astropy quantity) escape speed (km/s)
            r (astropy quantity) distance from the galactic center (kpc)
            a (astropy quantity) hernquist scaling (kpc) set at a default of 30 kpc
    output: M (astropy quantity) mass within r (Msun)
    '''
    vesckpcgyr = vesc.to(u.kpc/u.Gyr) #convert vesc into kpc/gyr
    M = vesckpcgyr**2 /2/Grav*(r+a) #calculates the mass
    return M

Vleo = 196*u.km/u.s #speed of Leo 1 satelite
r = 260*u.kpc
MLeo1 = MassHernVesc(Vleo,r)
print(f'Mass of Dark matter to keep LeoI bound: {MLeo1:.2e}')