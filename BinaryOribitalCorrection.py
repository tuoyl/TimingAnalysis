#!/usr/bin/env python

from __future__ import division
import numpy as np
from astropy.io import fits

__all__ = ['orbit_cor_bt', 'orbit_cor_deeter']

#correcting_method = ["Reggio", "Sanna", "DT"]
correcting_method = ["BT", "Deeter"]


def orbit_cor_bt(t, Porb, asini, e, omega, Tw, gamma):
    """
    use numerical method to solve Kepler equation and calculate delay
    BT model (Blandford & Teukolsky, 1976)

    Parameters
    -----------------
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

    gamma ; float
        the coefficient measures the combined effect of gravitational redshift and time dilation
        
    Returns 
    -------------
    new_t : array-like
        The time series after the binary correction

    """
    x = asini
    gamma = 0 

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
    factor1 = x*np.sin(omega)*(np.cos(E)-e) + (x*np.cos(omega)* ((1-e**2)**0.5) + gamma )*np.sin(E)
    #factor2
    factor2 = 1- (2*np.pi/Porb)* (x*np.cos(omega)*((1-e**2)**0.5)-x*np.sin(omega)*np.sin(E)) * (1-e*np.cos(E))**(-1)
    factor = factor1 * factor2
    new_t = t + factor #NOTE:pulsar proper Time = time + facotr
    print(factor)
    return new_t


def orbit_cor_deeter(time, Porb, asini, e, omega, Tnod):
    """
    Correct the photon arrival times to the photon emission time
    Deeter model (see, e.g., Deeter et al. 1981)

    Parameters
    -----------------
    time : array-like
        The time series before binary correction

    Porb : float
        The period of binary motion (in units of second)

    asini : float
        Projected semi-major orbital axis (in units of light-sec)

    e : float 
        The orbital eccentricity

    omega : float
        Longitude of periastron

    Tnod : float 
        The epoch of ascending node passage (in units of seconds, same time system with parameter t)

    Returns 
    -------------
    t_em : array-like
        The photon emission time

    """
    A = asini
    mean_anomaly = 2*np.pi*(time-Tnod)/Porb

    term1 = np.sin(mean_anomaly + omega) 
    term2 = (e/2)*np.sin(2*mean_anomaly + omega)
    term3 = (-3*e/2)*np.sin(omega)

    t_em = time - A * (term1 + term2 + term3)
    return t_em

if __name__ == "__main__":
    hdulist = fits.open("./testdata/P021100601401_he_screen.fits")
    time = hdulist[1].data.field("TDB")
