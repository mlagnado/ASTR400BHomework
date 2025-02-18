import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt

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


    def MassEnclosed(self,ptype,radius):

        COM = CenterOfMass(self.filename,ptype)

        COM_p = COM.COM_P(0.1)
        xR = self.x-COM_p[0]
        yR = self.y-COM_p[1]
        zR = self.z-COM_p[2]
        R = (xR**2 + yR**2 + zR**2)**0.5
        Masses = []
        for r in radius:
            indexR = np.where(R < r*u.kpc)
            m_tot = np.sum(self.m[indexR])
            Masses.append(m_tot*1e10*u.Msun)
        return Masses





if __name__ == '__main__' :
    MW = MassProfile("MW", 0)
    rs = np.arange(0.25,30.5,1.5)
    print(rs)
    MWmass = MW.MassEnclosed(1,rs)
    print(MWmass)
	##Cant plot units?
    plt.scatter(MWmass,rs)
