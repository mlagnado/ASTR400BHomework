import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt
from astropy import constants as const

class MassProfile:

    def __init__(self,galaxy,snap):
        '''
        This initializes the class getting the data from the file and creating the function name
        Inputs: 
            galaxy (string) should be the name of the galaxy ex. 'MW' for milky way
            snap (int) the snap number of the desired data
                       This i srelated to the timestep of the simulation.
        '''
        snap_num = '000'+str(snap) #Turns the snap number into a string
        snap_num = snap_num[-3:] #Makes sure the snap number is only 3 characters long
        self.filename = "%s_"%(galaxy) + snap_num + '.txt' #Connects the galaxy name and snap number to get the desired filename

        self.time, self.total, self.data = Read(self.filename) #Reads in data from the file
        self.x = self.data['x']*u.kpc #initializes x position data with units
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.m = self.data['m']
        self.gname = galaxy #collects the galaxy name for later use


    def MassEnclosed(self,ptype,radii):
        '''
        This function Calculates the mass enclosed by a desired radius and particle type
        Inputs:
            ptype (int) this relates to the particle type 1-halo, 2-disk, 3-bulge
            radii (list of floats) These are the radii where the mass enclosed will be calculated at.
        Outputs:
            A list of masses from each radius in units of solar masses
        '''
        COM = CenterOfMass(self.filename,ptype) #Calculates the center of mass
        COM_p = COM.COM_P(0.1) #Calculates the position coordinates of the COM
        indexP = np.where(self.data['type'] == ptype) #Finds the data of a specific particle type
        xR = self.x-COM_p[0] #Shifts the x position data by the COM coords
        yR = self.y-COM_p[1]
        zR = self.z-COM_p[2]
        xP = xR[indexP] #Filters the x position data by the particle type
        yP = yR[indexP]
        zP = zR[indexP]
        mP = self.m[indexP]
        R = (xP**2 + yP**2 + zP**2)**0.5 #Calculates the radius at each point
        Masses = np.zeros(shape=len(radii)) #Initializes a list of masses (Same length as list of radii)
        for i in range(len(radii)):
            indexR = np.where(R < radii[i]*u.kpc) #Finds points where the radius is less than the desired enclosed radius
            m_tot = np.sum(mP[indexR]) #Sums the mass inside the enclosed radius
            Masses[i] = m_tot #Adds the sum to the list
        Masses = Masses*1e10 #Adjusts the values for units
        return Masses*u.Msun #returns and adds units to the output

    def MassEnclosedTotal(self,radii):
        '''
        This function sums the mass of all particle types 
        Input:
            radii (list of floats) These are the radii where the mass enclosed will be calculated at
        Output:
            MTotal (list of astropy quantities) These are the total enclosed mass (units of Msun) at each input radius
        '''
        MHalo = self.MassEnclosed(1,radii) #calculates the mass enclosed by halo particle type
        MDisk = self.MassEnclosed(2,radii) #calculates the mass enclosed by disk particle type
        if self.gname == 'M33': #M33 does not have a bulge particle which we have to account for
            MTotal = MHalo + MDisk
        else: #If not M33 we calculate bulge enclosed mass
            MBulge = self.MassEnclosed(3,radii) #Calculates the mass enclosed by bulge particle type
            MTotal = MHalo + MDisk + MBulge
        return MTotal

    def HernquistMass(self,radii,a,Mhalo):
        '''
        This functions calculates the Hernquist Mass at each radius using the halo mass
        Inputs:
            radii (list of floats) These are the radii where the mass enclosed was calculated
            a (float) The hernquist scale radius 
            Mhalo (list of astropy quantities) The enclosed halo mass at each radius
        Output:
            halo_mass (list of floats) The calculated hernquist mass at each input radius
        '''
        numerator = Mhalo*(radii**2) #Calculates the numerator of the hernquist profile equation
        denominator = (a+radii)**2 #Calculates the denominator of the hernquist profile equation
        halo_mass = numerator/denominator #Calculates the hernquist mass profile
        return halo_mass

    def CircularVelocity(self,ptype,radii):
        '''
        This function calculates the circular velocity of each particle type at each radius
        Inputs:
            ptype (int) the desired particle type 1-halo,2-disk,3-bulge
            radii (list of floats) The radii at which the circular velocity of the mass enclosed will be calculated
        '''
        const.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun) #conversion of the G constant to desired units
        Mass = self.MassEnclosed(ptype,radii) #Calculates the mass enclosed at each radius
        circvel = np.sqrt(const.G*Mass/radii) #Calculates the circular velocity of the enclosed mass at each radius
        return circvel
    
    def CircularVelocityTotal(self,radii):
        '''
        This function calculates the total circular velocity of all particle types at each radius
        Input:
            radii (list of floats) The radius at which the circular velocity of all particle types will be calculated at
        Output:
            VTotal (list of floats) The calculated total circular velocity at each input radius
        '''
        const.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun) #Converts G constant to the desired units
        MTotal = self.MassEnclosedTotal(radii) #Calculates the total enclosed mass at each radius
        VTotal = np.sqrt(const.G*MTotal/radii) #Calculates the total circular velocity at each enclosed radius
        return VTotal

    def HernquistVCirc(self,radii,a,Mhalo):
        '''
        This function calculates the hernquist circular velocity
        Inputs:
            radius (list of floats) The radii at which the hernquist circular velocity will be calculated at
            a (float) The hernquist scale radius
            Mhalo (list of astropy quantities) The enclosed mass of the halo at each radius
        '''
        const.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun) #Converts the G constant into desired units
        Mass = self.HernquistMass(radii,a,Mhalo) #Calculates the hernquist mass at each radius
        VCirc = np.sqrt(const.G*Mass/radii) #Calculates the circular velocity of the hernquist mass at each radius
        return VCirc



if __name__ == '__main__' :
    '''
    Plots all the plots for the homework
    '''
    MW = MassProfile("MW", 0) #creates the mass profile for MW
    rs = np.arange(0.1,30.0,0.5) #List of radii used for all plots
    MWmassHalo = MW.MassEnclosed(1,rs) #Calculates the enclosed halo mass at each radius
    MWmassDisk = MW.MassEnclosed(2,rs) #Calculates the enclosed disk mass at each radius
    MWmassBulge = MW.MassEnclosed(3,rs) #Calculates the enclosed bulge mass at each radius
    MWmassTotal = MW.MassEnclosedTotal(rs) #Calculates the enclosed total mass at each radius
    a = 60 #Used in class for MW
    MWHernquist = MW.HernquistMass(rs,a,MWmassHalo) #Calculates the hernquist mass at each radius
    fig = plt.figure(figsize=(7,4))
    ax = plt.subplot(111)
    plt.plot(rs,MWmassHalo, label='Halo Mass Profile')
    plt.plot(rs,MWmassDisk, label='Disk Mass Profile')
    plt.plot(rs,MWmassBulge, label='Bulge Mass Profile')
    plt.plot(rs,MWmassTotal, label='Total Mass Profile')
    plt.plot(rs,MWHernquist, label=f'Hernquist Profile (a={a} kpc)')
    plt.yscale("log")
    plt.xlabel('Radius (kpc)')
    plt.title(f'Mass profile for MW')
    plt.ylabel(r'Log(Mass Enclosed ($M_{\odot}$))')
    plt.xlim(-1,32)
    plt.ylim(10**9,4*10**11)
    plt.legend()
    plt.show()
#    plt.savefig('MWMass.png')
    
    M31 = MassProfile("M31", 0) #Creates the mass profile for M31
    M31massHalo = M31.MassEnclosed(1,rs) #Calculates the enclosed halo mass at each radius
    M31massDisk = M31.MassEnclosed(2,rs)#Calculates the enclosed disk mass at each radius
    M31massBulge = M31.MassEnclosed(3,rs)#Calculates the enclosed bulge mass at each radius
    M31massTotal = M31.MassEnclosedTotal(rs) #Calculates the enclosed total mass at each radius
    a = 60 #Used because M31 is similar to MW
    M31Hernquist = M31.HernquistMass(rs,a,M31massHalo) #Calculates the hernquist mass at each radius
    ax = plt.subplot(111)
    plt.plot(rs,M31massHalo, label='Halo Mass Profile')
    plt.plot(rs,M31massDisk, label='Disk Mass Profile')
    plt.plot(rs,M31massBulge, label='Bulge Mass Profile')
    plt.plot(rs,M31massTotal, label='Total Mass Profile')
    plt.plot(rs,M31Hernquist, label=f'Hernquist Profile (a={a} kpc)')
    plt.yscale("log")
    plt.xlabel('Radius (kpc)')
    plt.title(f'Mass profile for M31')
    plt.ylabel(r'Log(Mass Enclosed ($M_{\odot}$))')
    plt.xlim(-1,32)
    plt.ylim(10**9,4*10**11)
    plt.legend()
    plt.show()
#    plt.savefig('M31Mass.png')

    M33 = MassProfile("M33", 0) #Creates the mas profile for M33
    M33massHalo = M33.MassEnclosed(1,rs) #Calculates the enclosed halo mass at each radius
    M33massDisk = M33.MassEnclosed(2,rs) #Calculates the enclosed disk mass at each radius
    #M33 does not have a Bulge
    M33massTotal = M33.MassEnclosedTotal(rs) #Calculates the enclosed total mass at each radius
    a = 25 #Used because M33 is smaller than M31 and MW
    M33Hernquist = M33.HernquistMass(rs,a,M33massHalo) #Calculates the hernquist mass at each radius
    ax = plt.subplot(111)
    plt.plot(rs,M33massHalo, label='Halo Mass Profile')
    plt.plot(rs,M33massDisk, label='Disk Mass Profile')
    plt.plot(rs,M33massTotal, label='Total Mass Profile')
    plt.plot(rs,M33Hernquist, label=f'Hernquist Profile (a={a} kpc)')
    plt.yscale("log")
    plt.xlabel('Radius (kpc)')
    plt.title(f'Mass profile for M33')
    plt.ylabel(r'Log(Mass Enclosed ($M_{\odot}$))')
    plt.xlim(-1,32)
    plt.ylim(10**9,4*10**11)
    plt.legend()
    plt.show()
#    plt.savefig('M33Mass.png')

    #plot 2
    a = 60 #Used in class for MW
    MWRotHalo = MW.CircularVelocity(1,rs) #Calculates the circular velocity for the enclosed halo mass
    MWRotDisk = MW.CircularVelocity(2,rs) #Calculates the circular velocity for the enclosed disk mass
    MWRotBulge = MW.CircularVelocity(3,rs) #Calculates the circular velocity for the enclosed bulge mass
    MWRotCurve=MW.CircularVelocityTotal(rs)  #Calculates the circular velocity for the enclosed total mass
    MWHernquistCirc = MW.HernquistVCirc(rs,a,MWmassHalo)  #Calculates the circular velocity for the hernquist mass
    ax = plt.subplot(111)
    plt.plot(rs,MWRotHalo, label='Halo Rotation Curve')
    plt.plot(rs,MWRotDisk, label='Disk Rotation Curve')
    plt.plot(rs,MWRotBulge, label='Bulge Rotation Curve')
    plt.plot(rs,MWRotCurve, label='Rotation Curve (Total)')
    plt.plot(rs,MWHernquistCirc, label=f'Hernquist Rotation Curve (a={a} kpc)')
    plt.title('Rotation Curve for MW')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Circular Velocity $(\frac{km}{s})$')
    plt.legend()
    plt.show()
#    plt.savefig('MWVel.png')

    a = 60
    M31RotHalo = M31.CircularVelocity(1,rs) #Calculates the circular velocity for the enclosed halo mass
    M31RotDisk = M31.CircularVelocity(2,rs) #Calculates the circular velocity for the enclosed disk mass
    M31RotBulge = M31.CircularVelocity(3,rs) #Calculates the circular velocity for the enclosed bulge mass
    M31RotCurve=M31.CircularVelocityTotal(rs) #Calculates the circular velocity for the enclosed total mass
    M31HernquistCirc = M31.HernquistVCirc(rs,a,M31massHalo) #Calculates the circular velocity for the hernquist mass
    ax = plt.subplot(111)
    plt.plot(rs,M31RotHalo, label='Halo Rotation Curve')
    plt.plot(rs,M31RotDisk, label='Disk Rotation Curve')
    plt.plot(rs,M31RotBulge, label='Bulge Rotation Curve')
    plt.plot(rs,M31RotCurve, label='Rotation Curve (Total)')
    plt.plot(rs,M31HernquistCirc, label=f'Hernquist Rotation Curve (a={a} kpc)')
    plt.title('Rotation Curve for M31')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Circular Velocity $(\frac{km}{s})$')
    plt.legend()
    plt.show()
#    plt.savefig('M31Vel.png')

    a = 30
    M33RotHalo = M33.CircularVelocity(1,rs) #Calculates the circular velocity for the enclosed halo mass
    M33RotDisk = M33.CircularVelocity(2,rs) #Calculates the circular velocity for the enclosed disk mass
    #M33 does not have a Bulge
    M33RotCurve=M33.CircularVelocityTotal(rs) #Calculates the circular velocity for the enclosed total mass
    M33HernquistCirc = M33.HernquistVCirc(rs,a,M31massHalo) #Calculates the circular velocity for the hernquist mass
    ax = plt.subplot(111)
    plt.plot(rs,M33RotHalo, label='Halo Rotation Curve')
    plt.plot(rs,M33RotDisk, label='Disk Rotation Curve')
    plt.plot(rs,M33RotCurve, label='Rotation Curve (Total)')
    plt.plot(rs,M33HernquistCirc, label=f'Hernquist Rotation Curve (a={a} kpc)')
    plt.title('Rotation Curve for M33')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Circular Velocity $(\frac{km}{s})$')
    plt.legend()
    plt.show()
#    plt.savefig('M33Vel.png')
    '''
    ##I could have shortened the plotting part of the code with a for loop
    but it was not that much extra to just copy paste it twice##

    I had to manually save these figures because they were coming back blank
    I tried to get the hernwuist profile to line up but could not find the error there
    '''