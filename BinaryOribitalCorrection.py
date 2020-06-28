#!/usr/bin/env python
"""
Definition of :class::class:`BinaryOrbitalCorrection`

:class::class:`BinaryOrbitalCorrection` is used to calculate the time of arrival
from the center of mass of the binary system.
"""

import numpy as np
from astropy.io import fits

__all__ = ["orbit_cor"]

#correcting_method = ["Reggio", "Sanna", "DT"]


def orbit_cor(t, Porb, asini, e, omega, Tw):
    """
    use numerical method to solve Kepler equation and calculate delay

    Parameters
    ____________
    t : array-like
        The time series without binary correction

    Porb : float
        The period of binary motion (in units of second)

    asini : float
        Projected semi-major orbital axis (in units of light-sec)

    e : float 
        The orbital eccentricity

    omega : float
        Longitude of periastron

    Tw : float 
        The epoch of periastron passage (in units of seconds, same time system with parameter t)
        
    Returns 
    ____________
    new_t : array-like
        The time series after the binary correction

    """
    x = asini

    if e == 0:
        E = 2*np.pi*(t-Tw)/Porb;
    else:
        E = np.array([])
        #solve Kepler Equation
        dE = 1e-3
        E_min = 2*np.pi*(t-Tw)/Porb - e
        E_max = 2*np.pi*(t-Tw)/Porb + e
        for i in range(len(E_min)): #TODO:the numerical method need to optimize!
            E_arr = np.arange(E_min[i], E_max[i], dE)

            equation_left = E_arr - e*np.sin(E_arr)
            equation_right= 2*np.pi*(t[i]-Tw)/Porb
            residual = np.abs(equation_left - equation_right)
            min_index = np.where(residual == min(residual))

            E = np.append(E, E_arr[min_index])
    
    #calculate time delay by orbit
    #factor1
    factor1 = x*np.sin(omega)*(np.cos(E)-e) + (x*np.cos(omega)* ((1-e**2)**0.5) )*np.sin(E)
    #factor2
    factor2 = 1- (2*np.pi/Porb)* (x*np.cos(omega)*((1-e**2)**0.5)-x*np.sin(omega)*np.sin(E)) * (1-e*np.cos(E))**(-1)
    factor = factor1 * factor2
    new_t = t + factor #NOTE:pulsar proper Time = time + facotr
    print(factor)
    return new_t


if __name__ == "__main__":
    hdulist = fits.open("./testdata/P021100601401_he_screen.fits")
    time = hdulist[1].data.field("TDB")
