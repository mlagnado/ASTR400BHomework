import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt
from astropy import constants as const

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
        indexP = np.where(self.data['type'] == ptype)
        xR = self.x-COM_p[0]
        yR = self.y-COM_p[1]
        zR = self.z-COM_p[2]
        xP = xR[indexP]
        yP = yR[indexP]
        zP = zR[indexP]
        mP = self.m[indexP]
        R = (xP**2 + yP**2 + zP**2)**0.5
        Masses = np.zeros(shape=len(radii))
        for i in range(len(radii)):
            indexR = np.where(R < radii[i]*u.kpc)
            m_tot = np.sum(mP[indexR])
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

    def HernquistMass(self,radii,a,Mhalo):
        numerator = Mhalo*(radii**2)
        denominator = (a+radii)**2
        halo_mass = numerator/denominator
        return halo_mass

    def CircularVelocity(self,ptype,radii):
        const.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        Mass = self.MassEnclosed(ptype,radii)
        circvel = np.sqrt(const.G*Mass/radii)
        return circvel
    
    def CircularVelocityTotal(self,radii):
        const.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        MTotal = self.MassEnclosedTotal(radii)
        VTotal = np.sqrt(const.G*MTotal/radii)
        return VTotal

    def HernquistVCirc(self,radius,a,Mhalo):
        const.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        Mass = self.HernquistMass(radius,a,Mhalo)
        VCirc = np.sqrt(const.G*Mass/radius)
        return VCirc



if __name__ == '__main__' :
    MW = MassProfile("MW", 0)
    rs = np.arange(0.1,30.0,0.5)
    MWmassHalo = MW.MassEnclosed(1,rs)
    MWmassDisk = MW.MassEnclosed(2,rs)
    MWmassBulge = MW.MassEnclosed(3,rs)
    MWmassTotal = MW.MassEnclosedTotal(rs)
    a = 40
    MWHernquist = MW.HernquistMass(rs,a,MWmassHalo)
    fig = plt.figure(figsize=(7,4))
    ax = plt.subplot(111)
    plt.plot(rs,MWmassHalo, label='Halo Mass Profile',alpha=0.5)
    plt.plot(rs,MWmassDisk, label='Disk Mass Profile',alpha=0.5)
    plt.plot(rs,MWmassBulge, label='Bulge Mass Profile',alpha=0.5)
    plt.plot(rs,MWmassTotal, label='Total Mass Profile', alpha=0.5)
    plt.plot(rs,MWHernquist, label='Hernquist Profile',alpha=0.5)
    plt.yscale("log")
    plt.xlabel('Radius (kpc)')
    plt.title(f'Mass profile for MW')
    plt.ylabel(r'Log(Mass Enclosed ($M_{\odot}$))')
    plt.xlim(-1,32)
    plt.ylim(10**9,4*10**11)
    plt.legend()
    plt.show()
    plt.savefig(f'MWMass.png')
    #####Create plots for M31 and M33##########

    #plot 2
    MWRotHalo = MW.CircularVelocity(1,rs)
    MWRotDisk = MW.CircularVelocity(2,rs)
    MWRotBulge = MW.CircularVelocity(3,rs)
    MWRotCurve=MW.CircularVelocityTotal(rs)
    MWHernquistCirc = MW.HernquistVCirc(rs,a,MWmassHalo)
    ax = plt.subplot(111)
    plt.plot(rs,MWRotHalo, label='Halo Rotation Curve', alpha=0.5)
    plt.plot(rs,MWRotDisk, label='Disk Rotation Curve', alpha=0.5)
    plt.plot(rs,MWRotBulge, label='Bulge Rotation Curve', alpha=0.5)
    plt.plot(rs,MWRotCurve, label='Rotation Curve (Total)', alpha=0.5)
    plt.plot(rs,MWHernquistCirc, label='Hernquist Rotation Curve', alpha=0.5)
    plt.title('')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Circular Velocity $(\frac{km}{s})$')
    plt.legend()
    plt.show()
    plt.savefig('MWVel.png')
    #####Create plots for M31 and M33##########
    #####Fix Hernquist Profile for both########