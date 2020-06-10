#!/usr/bin/env python
from __future__ import division
import numpy as np
from astropy.io import fits
from tqdm import tqdm
import argparse
import sys


def fsearch(filename, profilename, chisquarename, f0, fstep, frange, epoch, bin_profile, **kwargs):
    ''' NOTICE: this fsearch function is different with fsearch in function base in hxmtanalysis,
    which calculate best frequency without parameter file'''
    hdulist = fits.open(filename)
    try:
        time = hdulist[1].data.field("TDB")
    except:
        print("WARNING: the time series is read by column Time!")
        time = hdulist[1].data.field("TIME")
    
    MJDREFF = 0.0007660185
    MJDREFI = 55927

    #read parfile and parameters

    PEPOCH = epoch
    print("EPOCH = ", PEPOCH)
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    F0 = f0
    if kwargs["f1"] in kwargs:
        f1 = kwargs["f1"]
    else:
        f1 = 0
    if kwargs["f2"] in kwargs:
        f2 = kwargs["f2"]
    else:
        f2 = 0
    if kwargs["f3"] in kwargs:
        f3 = kwargs["f3"]
    else:
        f3 = 0


    tf = np.array([])
    fre = np.array([])
    fre_err = np.array([])

    
    data = time
    if len(data)==0:
        raise IOError("Error: Data is empty")
    #t0 = min(data)
    t0 = pepoch
    #T0 = t0/86400.0 + MJDREFF + MJDREFI
    #dt = t0 - pepoch 
    #f0 = F0 + F1*dt + (1/2)*F2*(dt**2) + (1/6)*F3*(dt**3) + (1/24)*F4*(dt**4)
    #f1 = F1 + F2*dt + (1/2)*F3*(dt**2) + (1/6)*F4*(dt**4)
    #f2 = F2 + F3*dt + (1/2)*F4*(dt**2)
    #f3 = F3 + F4*dt
    #f4 = F4

    f = np.arange(f0-frange,f0+frange,fstep)
    chi_square = np.array([])
    N = len(data)
    b = N/bin_profile
    #for f1 in np.arange(-3.72e-10,-3.69e-10,0.5e-12):
    for j in tqdm(range(0,len(f))):
        phi_tmp = (data-t0)*f[j] + (1.0/2)*((data-t0)**2)*f1 + (1.0/6.0)*((data-t0)**3)*f2 + (1.0/24)*((data-t0)**4)*f3
        phi_tmp -= np.floor(phi_tmp)
        p_num = np.histogram(phi_tmp,bin_profile)[0]
        #chi_square[j] = np.std(p_num)**2/np.mean(p_num)
        chi_square = np.append(chi_square,(np.std(p_num)**2/np.mean(p_num)))

    print '\n'
    fbest = f[np.where(chi_square==max(chi_square))][0]

    #save chisquare
    with open(chisquarename, 'w') as fout:
        for i in range(len(f)):
            fout.write("%.12f %.6f \n"%(f[i], chi_square[i]))
    
    #profiles
    print(data, t0, fbest)
    phi_tmp = (data-t0)*fbest + (1.0/2)*((data-t0)**2)*f1 + (1.0/6.0)*((data-t0)**3)*f2 + (1.0/24)*((data-t0)**4)*f3
    phi_tmp -= np.floor(phi_tmp)
    profile_best, profile_x = np.histogram(phi_tmp,bin_profile)
    profile_x = profile_x[:-1]
    with open(profilename, 'w')as fout:
        for i in range(len(profile_x)):
            fout.write("%f %f %f\n"%(profile_x[i], profile_best[i], np.sqrt(profile_best[i])))


if __name__ == "__main__" :
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Example: python he_pipeline.py -i hxmt_filename -p outprofile.dat -c outchisquare.dat')
    parser.add_argument("-i","--input",help="input filename of HXMT screen file")
    parser.add_argument("-p","--profile",help="out profile name")
    parser.add_argument("-c","--chisquare",help="chisquare distribution file",type=str)
    parser.add_argument("-f0","--freqency", help="f0 for searching frequency", type=float)
    parser.add_argument("-f1","--freqderive", help="f1 for searching frequency", type=float)
    parser.add_argument("-f2","--freqsecderive", help="f2 for searching frequency", type=float)
    parser.add_argument("-f3","--freqthirdderive", help="f3 for searching frequency", type=float)
    parser.add_argument("-fstep",help="frequency intervals for frequency search", type=float)
    parser.add_argument("-frange", help="frequency range for searching", type=float)
    parser.add_argument("-epoch", help="epoch time for fsearch (MJD)", type=float)
    parser.add_argument("-bins", help="profile bin numbers", type=int)
    args = parser.parse_args()
    filename = args.input
    outprofile = args.profile
    outchisq   = args.chisquare
    f0 = args.freqency
    fstep = args.fstep
    frange = args.frange
    epoch  = args.epoch
    binprofile = args.bins
    if args.freqderive:
        f1 = args.freqderive
    else:f1 = 0
    if args.freqsecderive:
        f2 = args.freqsecderive
    else:f2 = 0
    if args.freqthirdderive:
        f3 = args.freqthirdderive
    else:f3 = 0

    fsearch(filename, outprofile, outchisq, f0, fstep, frange, epoch, binprofile, f1=f1, f2=f2, f3=f3)
