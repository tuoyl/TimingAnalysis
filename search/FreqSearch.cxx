#include <iostream>
#include <algorithm>
#include <string.h>
#include "fitsio.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include "/Users/tuoyouli/cpplibrary/MyFun.h"
#include "../include/matplotlibcpp.h"
namespace plt = matplotlibcpp;

//prototype
void PrintError(int status);
char* getCmdOption(char** begin, char** end, const std::string& option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
//

/* global variables */
/* TODO: read parameters from parfile */ 
double MinFreq = 11.0846499395;
double MaxFreq = 11.3046499395;
double FreqStep= 1e-6;
int MJDREFI = 51910;
double MJDREFF = 0.00074287037037037;
double PEpoch  = (51559.319 - MJDREFI -MJDREFF)*86400;
double F1 = 0;
double F2 = 0;
int bin_cs = 20;

int main(int argc, char* argv[])
{
    /* READ parse arguments */
    if(cmdOptionExists(argv, argv+argc, "-h"))
    {
        std::cout << "########## HELP ############" << std::endl;
    }
//    char* InputEventFile = argv[1];
    char* InputEventFile = getCmdOption(argv, argv + argc, "-f");
    if (InputEventFile)
    {
        std::cout << "############################" << std::endl;
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
        // NOTE: if "TDB" exist read "TDB", if not use "Time".
        
        if (fits_get_colnum(fptr, CASEINSEN, "Time", &ColnumTime, &status)) PrintError(status);
        if (ffgky(fptr, TSTRING, "INSTRUME", &InstrumentName, NULL, &status)) PrintError(status);
        printf("Instrument is '%s'. \n", InstrumentName);
    
    
    
        std::cout << "Column number is " << ColnumTime << std::endl;
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
        plt::plot(FreqArray, Chisquare, "o-");
        plt::show();

        std::cout << std::endl << "FINISH Frequency Search" << std::endl;
    
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





