#include <iostream>
#include <algorithm>
#include <string.h>
#include "fitsio.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "/Users/tuoyouli/cpplibrary/MyFun.h"
#include "../include/matplotlibcpp.h"
namespace plt = matplotlibcpp;

//prototype
void PrintError(int status);
char* getCmdOption(char** begin, char** end, const std::string& option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
void ReadParFile(std::string ParfileName, double& PEPOCH, double& F0, double& FreqRange, double& FreqStep, int& BinNum, std::string& TimeColumnName);
//

///* global variables */
double PEPOCH;
double PEpoch;
int MJDREFI;
double MJDREFF;
double F0 = 0;
double F1 = 0;
double F2 = 0;
double FreqStep;
double FreqRange;
double MinFreq;
double MaxFreq;
int bin_cs;
std::string TimeColumnName;

int main(int argc, char* argv[])
{
    /* READ parse arguments */
    if(cmdOptionExists(argv, argv+argc, "-h"))
    {
        std::cout << "########## HELP ############" << std::endl << std::endl;
        std::cout << "USAGE: " << std::endl;
        std::cout << "        FreqSearch -f EventFile.Fits -p ParameterFile.par -o SaveProfile.FITS" << std::endl << std::endl;
        std::cout << "########## END HELP ############" << std::endl;
    }
    char* InputEventFile = getCmdOption(argv, argv + argc, "-f");
    char* ParameterFile  = getCmdOption(argv, argv + argc, "-p");
    char* OutputProfileFile = getCmdOption(argv, argv + argc, "-o");

    /* READ pulsar parameters */
    if (ParameterFile)
    {
        ReadParFile(ParameterFile, PEPOCH, F0, FreqRange, FreqStep, bin_cs, TimeColumnName);
        if (F0-FreqRange < 0) 
        {
            MinFreq = 0;
        }
        else
        {
            MinFreq = F0 - FreqRange;
        }
        MaxFreq = F0 + FreqRange;
        printf("F0: %f\n", F0);
        printf("FreqRange: %f\n", FreqRange);
    }

    if (InputEventFile)
    {
        std::cout << "########## Searching ############" << std::endl;
        std::cout << InputEventFile << std::endl;
    
        /* Read the Time of fits file */
        fitsfile* fptr;
        int status = 0;
        int datatype, anynull;
        double doublenull = 0;
        int LengthArrayTime;
        // get length of array
        if (fits_open_file(&fptr, InputEventFile, READONLY, &status)) PrintError(status); 
        if (ffmahd(fptr, 2, &datatype, &status)) PrintError(status);
        if (ffgky(fptr, TINT, "NAXIS2", &LengthArrayTime, NULL, &status)) PrintError(status);
    
        double* TimeArray;
        TimeArray = new double[LengthArrayTime];
    
        int ColnumTime;
        char InstrumentName[10]; 
        char TimeColumnCharName[TimeColumnName.length()];
        for (int i=0; i < sizeof(TimeColumnName); i++) 
        {
            TimeColumnCharName[i] = TimeColumnName[i];
        }
        
        if (fits_get_colnum(fptr, CASEINSEN, TimeColumnCharName, &ColnumTime, &status)) PrintError(status);
        if (ffgky(fptr, TSTRING, "INSTRUME", &InstrumentName, NULL, &status)) PrintError(status);
        printf("Instrument is '%s'. \n", InstrumentName);

        // GET THE MJD reference of Instrument
        if (ffgky(fptr, TINT, "MJDREFI", &MJDREFI, NULL, &status)) PrintError(status);
        std::cout << "MJDREFI = " << MJDREFI << std::endl;
        if (ffgky(fptr, TDOUBLE, "MJDREFF", &MJDREFF, NULL, &status)) PrintError(status);
        std::cout << "MJDREFF = " << MJDREFF << std::endl;
        double PEpoch = (PEPOCH - MJDREFI -MJDREFF)*86400; // PEpoch is MET time and PEPOCH is MJD time
    
        printf("Column number is %d ('%s') \n", ColnumTime, TimeColumnCharName);
        if (fits_read_col(fptr, TDOUBLE, ColnumTime, 1, 1, LengthArrayTime, &doublenull, TimeArray, &anynull, &status)) PrintError(status);
        // convert array to vector
        std::vector<double> TimeVector;
        for (int i=0; i<LengthArrayTime; i++)
        {
            TimeVector.push_back(TimeArray[i]);
        }
    
        /* MAIN algorithm  */
        std::vector<double> FreqArray;
        FreqArray = MyFun::Arange(MinFreq, MaxFreq, FreqStep);
        std::vector<double> Chisquare;
        std::vector<double> Phi;
        long double percentage;
        for (int i=0; i<FreqArray.size(); i++)
        {
            percentage = i/(double) FreqArray.size() * 100;
            std::cout.precision(4);
            std::cout << std::fixed << "\r" << percentage << "%" << std::flush;
            for (int j=0; j<TimeVector.size(); j++)
            {
                double Phitmp = (TimeVector[j] - PEpoch)*FreqArray[i] +
                        0.5*std::pow((TimeVector[j] - PEpoch), 2)*F1 +
                        (1.0/6.0)*std::pow((TimeVector[j] - PEpoch), 3)*F2;
                double Phifrac, Phiint;
                Phifrac = std::modf(Phitmp, &Phiint);
                Phi.push_back(Phifrac);
            }
            // NOTE: Calculate Chi-Square
            std::vector<int> ProfileCounts;
            ProfileCounts = MyFun::Hist1D(Phi, bin_cs);
            double chi2_i=0;
            double Ei = TimeVector.size()/(double) bin_cs;
            for (int k=0; k<ProfileCounts.size(); k++)
            {
                //Calculate chisquare test, profile counts is not actual profile counts here
                chi2_i += std::pow((ProfileCounts[k] - Ei),2)/Ei;
            }
            Chisquare.push_back(chi2_i);
    
            // clear Vector
            Phi.clear();
            ProfileCounts.clear();
        }
        for (int i; i<Chisquare.size(); i++)
        {std::cout << Chisquare[i] << " " << std::endl;}
        plt::plot(FreqArray, Chisquare, "o-");
        plt::show();
    }

    return 0;
}

void PrintError(int status)
{
    if (status)
    {
        fits_report_error(stderr, status);
        exit(status);
    }
    return;
}

char* getCmdOption(char** begin, char** end, const std::string& option)
{
    char** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}





