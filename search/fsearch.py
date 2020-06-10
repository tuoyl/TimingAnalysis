#!/usr/bin/env python
from __future__ import division
import numpy as np
from astropy.io import fits
from tqdm import tqdm
import argparse
import sys
import matplotlib.pyplot as plt


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
    if "f1" in kwargs:
        f1 = kwargs["f1"]
        if "f1step" in kwargs:
            f1step = kwargs["f1step"]
            f1range = kwargs["f1range"]
    else:
        f1 = 0
        f1step = 0 
        frange = 0
    if "f2" in kwargs:
        f2 = kwargs["f2"]
    else:
        f2 = 0
    if "f3" in kwargs:
        f3 = kwargs["f3"]
    else:
        f3 = 0


    tf = np.array([])
    fre = np.array([])
    fre_err = np.array([])

    
    data = time
    if len(data)==0:
        raise IOError("Error: Data is empty")
    t0 = pepoch

    f = np.arange(f0-frange,f0+frange,fstep)
    if f1 == 0:
        f1search = f1
    elif "f1step" in kwargs:
        f1search = np.arange(f1-f1range, f1+f1range, f1step)
    else:
        f1search = f1
        
    chi_square = np.array([])
    N = len(data)
    b = N/bin_profile


    if "f1range" in kwargs:
        print("F1 search space: ", f1search)
        print("F search space: ", f)
        chi_square = np.zeros((len(f1search), len(f)))
        for k in tqdm(range(0,len(f1search))):
            for j in range(len(f)):
                phi_tmp = (data-t0)*f[j] + (1.0/2)*((data-t0)**2)*f1search[k] + (1.0/6.0)*((data-t0)**3)*f2 + (1.0/24)*((data-t0)**4)*f3
                phi_tmp -= np.floor(phi_tmp)
                p_num = np.histogram(phi_tmp,bin_profile)[0]
                chi_square[k][j] = np.std(p_num)**2/np.mean(p_num)

    
        print '\n'
        maxchi2_index = np.where(chi_square == np.ndarray.max(chi_square))
        fbest = f[maxchi2_index[1][0]]
        f1best = f1search[maxchi2_index[0][0]]
        print(fbest, f1best)
    
        #save chisquare
        with open(chisquarename, 'w') as fout:
            for i in range(len(f1search)):
                for j in range(len(f)):
                    fout.write("%.12f %.12f %.6f \n"%(f[j], f1search[i], chi_square[i][j]))
        
        #profiles
        print(data, t0, fbest)
        phi_tmp = (data-t0)*fbest + (1.0/2)*((data-t0)**2)*f1best + (1.0/6.0)*((data-t0)**3)*f2 + (1.0/24)*((data-t0)**4)*f3
        phi_tmp -= np.floor(phi_tmp)
        profile_best, profile_x = np.histogram(phi_tmp,bin_profile)
        profile_x = profile_x[:-1]
        with open(profilename, 'w')as fout:
            for i in range(len(profile_x)):
                fout.write("%f %f %f\n"%(profile_x[i], profile_best[i], np.sqrt(profile_best[i])))

    # ONLY SEARCH F0 parameter space
    else:
        print("F search space: ", f)
        chi_square = np.zeros(len(f))
        for j in tqdm(range(0,len(f))):
            phi_tmp = (data-t0)*f[j] + (1.0/2)*((data-t0)**2)*f1 + (1.0/6.0)*((data-t0)**3)*f2 + (1.0/24)*((data-t0)**4)*f3
            phi_tmp -= np.floor(phi_tmp)
            p_num = np.histogram(phi_tmp,bin_profile)[0]
            chi_square[j] = np.std(p_num)**2/np.mean(p_num)
    
        ########
    
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
    parser.add_argument("-f1step",help="frequency derivative intervals for frequency search", type=float)
    parser.add_argument("-f1range", help="frequency derivative range for searching", type=float)
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
        freqderive = args.freqderive
    else:freqderive = 0
    if args.freqsecderive:
        freqsecderive = args.freqsecderive
    else:freqsecderive = 0
    if args.freqthirdderive:
        freqthirdderive = args.freqthirdderive
    else:freqthirdderive = 0

    if args.f1range:
        print("DO THIS")
        f1range = args.f1range
        f1step  = args.f1step
        fsearch(filename, outprofile, outchisq, f0, fstep, frange, epoch, binprofile, f1=freqderive, f2=freqsecderive, f3=freqthirdderive, f1step=f1step, f1range=f1range)
    else:
        fsearch(filename, outprofile, outchisq, f0, fstep, frange, epoch, binprofile, f1=freqderive, f2=freqsecderive, f3=freqthirdderive)
