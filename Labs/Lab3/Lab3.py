
#In Class Lab 3 Template
# G Besla ASTR 400B

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


# The Figure illustrates the color magnitude diagram (CMD) for the Carina Dwarf along with the interpreted 
# star formation history from isochrone fitting to the CMD.
# The image is from Tolstoy+2009 ARA&A 47 review paper about dwarf galaxies
# 
# ![Iso](./Lab3_Isochrones.png)
# opened outside of python

# # This Lab:
# 
# Modify the template file of your choice to plot isochrones that correspond to the inferred star formation episodes (right panel of Figure 1) to recreate the dominant features of the CMD of Carina (left panel of Figure 1). 



# Some Notes about the Isochrone Data
# DATA From   http://stellar.dartmouth.edu/models/isolf_new.html
# files have been modified from download.  ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading
# NOTE SETTINGS USED:  Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# These could all be changed and it would generate a different isochrone




# Filename for data with Isochrone fit for 1 Gyr
# These files are located in the folder IsochroneData
#filename1="./IsochroneData/Isochrone1.txt"
#filename2="./IsochroneData/Isochrone2.txt"
#filename3="./IsochroneData/Isochrone3.txt"
#filename4="./IsochroneData/Isochrone4.txt"
#filename5="./IsochroneData/Isochrone5.txt"
#filename6="./IsochroneData/Isochrone6.txt"
#filename7="./IsochroneData/Isochrone7.txt"
#filename8="./IsochroneData/Isochrone8.txt"
#filename9="./IsochroneData/Isochrone9.txt"
#filename10="./IsochroneData/Isochrone10.txt"
#filename11="./IsochroneData/Isochrone11.txt"
#filename12="./IsochroneData/Isochrone12.txt"
#filename13="./IsochroneData/Isochrone13.txt"
filename13_5="./IsochroneData/Isochrone13_5.txt"
#Can also do it as a for loop

# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

'''
I read all the data and plot it later to have one for loop
'''

# Read in data for an isochrone corresponding to 1 Gyr
#data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=8)
#names=true means we include column headers
#data2 = np.genfromtxt(filename2,dtype=None,names=True,skip_header=8)
#data3 = np.genfromtxt(filename3,dtype=None,names=True,skip_header=8)
#data4 = np.genfromtxt(filename4,dtype=None,names=True,skip_header=8)
#data5 = np.genfromtxt(filename5,dtype=None,names=True,skip_header=8)
#data6 = np.genfromtxt(filename6,dtype=None,names=True,skip_header=8)
#data7 = np.genfromtxt(filename7,dtype=None,names=True,skip_header=8)
#data8 = np.genfromtxt(filename8,dtype=None,names=True,skip_header=8)
#data9 = np.genfromtxt(filename9,dtype=None,names=True,skip_header=8)
#data10 = np.genfromtxt(filename10,dtype=None,names=True,skip_header=8)
#data11 = np.genfromtxt(filename11,dtype=None,names=True,skip_header=8)
#data12 = np.genfromtxt(filename12,dtype=None,names=True,skip_header=8)
#data13 = np.genfromtxt(filename13,dtype=None,names=True,skip_header=8)
data13_5 = np.genfromtxt(filename13_5,dtype=None,names=True,skip_header=8)


# Plot Isochrones 
# For Carina

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot Isochrones
for i in range(1,14):
    filename = f"./IsochroneData/Isochrone{i}.txt"
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=8)
    plt.plot(data['B']-data['R'], data['R'], linewidth=5, label=f"{i} Gyr")
# Isochrone for 1 Gyr
# Plotting Color vs. Difference in Color 
#plt.plot(data1['B']-data1['R'], data1['R'], color='blue', linewidth=5, label='1 Gyr')
###EDIT Here, following the same format as the line above 
#plt.plot(data2['B']-data2['R'], data2['R'], linewidth=5, label='2 Gyr')
#plt.plot(data3['B']-data3['R'], data3['R'], linewidth=5, label='3 Gyr')
#plt.plot(data4['B']-data4['R'], data4['R'], linewidth=5, label='4 Gyr')
#plt.plot(data5['B']-data5['R'], data5['R'], linewidth=5, label='5 Gyr')
#plt.plot(data6['B']-data6['R'], data6['R'], linewidth=5, label='6 Gyr')
#plt.plot(data7['B']-data7['R'], data7['R'], linewidth=5, label='7 Gyr')
#plt.plot(data8['B']-data8['R'], data8['R'], linewidth=5, label='8 Gyr')
#plt.plot(data9['B']-data9['R'], data9['R'], linewidth=5, label='9 Gyr')
#plt.plot(data10['B']-data10['R'], data10['R'], linewidth=5, label='10 Gyr')
#plt.plot(data11['B']-data11['R'], data11['R'], linewidth=5, label='11 Gyr')
#plt.plot(data12['B']-data12['R'], data12['R'], linewidth=5, label='12 Gyr')
#plt.plot(data13['B']-data13['R'], data13['R'], linewidth=5, label='13 Gyr')
plt.plot(data13_5['B']-data13_5['R'], data13_5['R'], linewidth=5, label='13.5 Gyr')

# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel('M$_R$', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.6, 0.15, 'CMD for Carina dSph', fontsize=22)

plt.savefig('IsochroneCarina.png')


print('## Question 2 ##')
# # Q2
# 
# Could there be younger ages than suggested in the Tolstoy plot?
# Try adding younger isochrones to the above plot.
# 
print('It would be possible to have ages younger than 1 Gyr')
#It would be possible to have ages younger than 1 Gyr

print('## Question 3 ##')
# # Q3
# 
# What do you think might cause the bursts of star formation?
# 
print('Bursts of star formation might be caused by collisions')
#Bursts of star formation might be caused by collisions
#If gas is too hot gravity cant pull everything together to cause star formation
#~10K or less