import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt
from astropy.constants import G

class MassProfile:

    def __init__(self,galaxy,snap):
        snap_num = '000'+str(snap)
        snap_num = snap_num[-3:]
        self.filename = "%s_"%(galaxy) + snap_num + '.txt'

        self.time, self.total, self.data = Read(self.filename)
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.m = self.data['m']
        self.gname = galaxy


    def MassEnclosed(self,ptype,radii):

        COM = CenterOfMass(self.filename,ptype)

        COM_p = COM.COM_P(0.1)
        xR = self.x-COM_p[0]
        yR = self.y-COM_p[1]
        zR = self.z-COM_p[2]
        R = (xR**2 + yR**2 + zR**2)**0.5
        Masses = np.zeros(shape=len(radii))
        for i in range(len(radii)):
            indexR = np.where(R < radii[i]*u.kpc)
            m_tot = np.sum(self.m[indexR])
            Masses[i] = m_tot
        Masses = Masses*1e10
        return Masses*u.Msun

    def MassEnclosedTotal(self,radii):
        MHalo = self.MassEnclosed(1,radii)
        MDisk = self.MassEnclosed(2,radii)
        if self.gname == 'M33':
            MTotal = MHalo + MDisk
        else:
            MBulge = self.MassEnclosed(3,radii)
            MTotal = MHalo + MDisk + MBulge
        return MTotal

    def HernquistMass(self,radius,a,Mhalo):
        numerator = Mhalo*radius**2
        denominator = (a+radius)**2
        halo_mass = numerator/denominator
        return halo_mass

    def CicularVelocity(self,ptype,radii):
        circvel = np.zeros(shape=len(radii))
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        for i in range(len(radii)):
            Mass = self.MassEnclosed(ptype,radii[i])
            circvel[i] = np.sqrt(G*Mass/radii[i])
        return circvel
    
    def CircularVelocityTotal(self,radii):
        VTotal = np.zeros(shape=len(radii))
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        for i in range(len(radii)):
            MTotal = self.MassEnclosedTotal(radii[i])
            VTotal[i] = np.sqrt(G*MTotal/radii[i])
        return VTotal

    def HernquistVCirc(self,radius,a,Mhalo):
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        Mass = self.HernquistMass(radius,a,Mhalo)
        VCirc = np.sqrt(G*Mass/radius)
        return VCirc



if __name__ == '__main__' :
    MW = MassProfile("MW", 0)
    rs = np.arange(0.25,30.5,.5)
    MWmass = MW.MassEnclosed(1,rs) 
    fig = plt.figure(figsize=(11,8))
    ax = plt.subplot(111)
    plt.scatter(rs,MWmass)
    plt.yscale("log")
    plt.xlabel('Radius (kpc)')
    plt.title(f'Mass profile for MW Halo')
    plt.ylabel(r'Log(Mass Enclosed ($M_{\odot}$))')
    plt.xlim(-1,32)
    plt.ylim(10**9,4*10**11)
    plt.savefig('P2plot.png')