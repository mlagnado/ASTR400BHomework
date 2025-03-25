
# # Lab 8 : Star Formation 



import numpy as np
from astropy import units as u
from astropy import constants as const

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#get_ipython().run_line_magic('matplotlib', 'inline')

print('## Part A ##')
# # Part A
# 
# Create a function that returns the SFR for a given luminosity (NUV, FUV, TIR, Halpha)
# 
# $Log( {\rm SFR} (M_\odot/year)) = Log(Lx (erg/s)) - Log(Cx)$ 
# 
# Including corrections for dust absorption 
# 
# Kennicutt & Evans 2012 ARA&A Equation 12 and Table 1, 2

#define our function for Star Formation Rate
def StarFormationRate(L, Type, TIR=0):
    '''
    Fxn that computes the SFR of a galaxy following the Kennicutt & Evans 2012 paper (eq 12) {ARA&A}
    Input param:
        L: float, luminosity of a galaxy in {ergs/sec}
        Type: string, tells the wavelength {'FUV','NUV','TIR','Halpha'}
        TIR: float, total infrared luminosity {ergs/sec} default value is 0
    Output:
        SFR: float, log of the star formation rate {Msun/yr}
    '''
    if Type == 'FUV':
        logCx = 43.35 #calibration from LFUV to SFR from the table 1 in Kennicutt & Evans paper
        TIRc = 0.46 #correction for dust absorption from table2 in K&E 2012 paper

    elif Type == 'NUV':
        logCx = 43.17 #calibration from LNUV to SFR from the table 1 in Kennicutt & Evans paper
        TIRc = 0.27 #correction for dust absorption from table2 in K&E 2012 paper

    elif Type == 'Halpha':
        logCx = 41.27 #calibration from LHalpha to SFR from the table 1 in Kennicutt & Evans paper
        TIRc = 0.0024 #correction for dust absorption from table2 in K&E 2012 paper

    elif Type == 'TIR':
        logCx = 43.41 #calibration from LFUV to SFR from the table 1 in Kennicutt & Evans paper
        TIRc = 0.0 #no correction

    else:
        print('Invalid Input:\nMissing wavelength input!!!\nWavelength input should be one of these {\'FUV\',\'NUV\',\'TIR\',\'Halpha\'}')

    #Correct the luminosity for dust using TIR
    Lcorr = L + TIRc*TIR

    #Star Formation rate
    SFR = np.log10(Lcorr)-logCx

    return SFR

#Testing if breaking the function works :)
#StarFormationRate(0.3,'a')

# Let's try to reproduce SFRs derived for the WLM Dwarf Irregular Galaxy using UV luminosities measured with Galex. 
# 
# Compare results to Table 1 from Lee et al. 2009 (who used the older Kennicutt 98 methods)
# https://ui.adsabs.harvard.edu/abs/2009ApJ...706..599L/abstract
# 
# We will use galaxy properties from NED (Photometry and SED):
# https://ned.ipac.caltech.edu/



# First need the Luminosity of the Sun in the right units
#astropy gives L_sun in watts but we need erg/sec
LsunErgS = const.L_sun.to(u.erg/u.s)
print(LsunErgS)

#  WLM Dwarf Irregular Galaxy
NUV_WLM = 1.71e7*LsunErgS #Value obtained from NED database for WLM galaxy
TIR_WLM = 3.21e5*LsunErgS + 2.49e6*LsunErgS + 2.48e6*LsunErgS #midIR+farIR+nearIR

#test SFR functions
testSFR = StarFormationRate(NUV_WLM.value, 'NUV', TIR_WLM.value)
print(testSFR)


print('## Part B ##')
# # Part B Star formation main sequence
# 
# 1) Write a function that returns the average SFR of a galaxy at a given redshift, given its stellar mass
# 
# 2) What is the average SFR of a MW mass galaxy today? at z=1?
# 
# 3) Plot the SFR main sequence for a few different redshifts from 1e9 to 1e12 Msun.
# 
# 
# From Whitaker 2012:
# 
# log(SFR) = $\alpha(z)({\rm log}M_\ast - 10.5) + \beta(z)$
# 
# $\alpha(z) = 0.7 - 0.13z$
# 
# $\beta(z) = 0.38 + 1.14z - 0.19z^2$

# # Step 1

#write a function that returns the average SFR of a galaxy at a given redshift given its stellar mass
def SFRMainSequence(Mstar, z):
    '''
    This is a function that computes the avg SFR of a galaxy
    Inputs:
        Mstar: float, stellar mass of the galaxy {Msun}
        z: float, redshift
    Outputs:
        SFR: float, log of the Star Formatino Rate {Msun} 
    '''
    alpha=0.7-0.13*z
    beta = 0.38 + 1.14*z -0.19*z**2
    SFR = alpha*(np.log10(Mstar) -10.5)+beta
    return SFR

# # Step 2



# MW at z=0




# MW at z = 1


# # Step 3



# create an array of stellar masses





fig = plt.figure(figsize=(8,8), dpi=500)
ax = plt.subplot(111)

# add log log plots


# Add axis labels
plt.xlabel('Log(Mstar (M$_\odot$))', fontsize=12)
plt.ylabel('Log(SFR (M$_\odot$/year))', fontsize=12)


#adjust tick label font size
label_size = 12
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')


# Save file
#plt.savefig('Lab8_SFR_MainSequence.png')

print('## Part C ##')
# # Part C  Starbursts
# 
# Use your `StarFormationRate` code to determine the typical star formation rates for the following systems with the listed Total Infrared Luminosities (TIR): 
# 
# Normal Galaxies: $10^{10}$ L$_\odot$
# 
# LIRG: $10^{11}$ L$_\odot$
# 
# ULIRG: $10^{12} $ L$_\odot$
# 
# HLIRG: $10^{13} $ L$_\odot$



# normal galaxies 




# LIRGs  




# ULIRGs




# HLIRGs






