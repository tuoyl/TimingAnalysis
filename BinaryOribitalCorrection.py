#!/usr/bin/env python
"""
Definition of :class::class:`BinaryOrbitalCorrection`

:class::class:`BinaryOrbitalCorrection` is used to calculate the time of arrival
from the center of mass of the binary system.
"""

import numpy as np
from astropy.io import fits

__all__ = ["BinaryOrbitalCorrection"]

correcting_method = ["Reggio", "Sanna", "DT"]

class BinaryOrbitalCorrection(object):

    def __init__(self, time, axsini=None, T_halfpi=None,
            P_orb=None, omega=None, e=None):
        self.time  = time        # time array in units of second
        self.axsini = axsini     # 
        self.T_halfpi = T_halfpi
        self.P_orb = P_orb
        self.omega = omega
        self.e = e


if __name__ == "__main__":
    hdulist = fits.open("./testdata/P021100601401_he_screen.fits")
    time = hdulist[1].data.field("TDB")
    toa = BinaryOrbitalCorrection(time)

    toa.P_orb =    22
    toa.axsini=    104.343 
    toa.P_orb=     22.534777 
    toa.omega=     220.10 
    toa.e    =     0.0150 
    toa.T_halfpi=  55927.26871 

    print(toa.P_orb)
    print(toa.e)
    

