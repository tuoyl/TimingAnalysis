#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include "/Users/tuoyouli/cpplibrary/MyFun.h"

//prototype
void PrintError(int status);
//

/* global variables */
/* TODO: read parameters from parfile */ 
double MinFreq = 26;
double MaxFreq = 28;
double FreqStep= 1e-2;
int MJDREFI = 51910;
double MJDREFF = 0.00074287037037037;
double PEpoch  = (51559.319 - MJDREFI -MJDREFF)*86400;
double F1 = 0;
double F2 = 0;
int bin_cs = 20;

int main(int argc, char* argv[])
{
    char* InputEventFile = argv[1];
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
        for (int i=0; i<TimeVector.size(); i++)
        {
            double Phitmp = (TimeVector[i] - PEpoch)*FreqArray[i] +
                    0.5*std::pow((TimeVector[i] - PEpoch), 2)*F1 +
                    (1.0/6.0)*std::pow((TimeVector[i] - PEpoch), 3)*F2;
            double Phifrac, Phiint;
            Phifrac = std::modf(Phitmp, &Phiint);
            Phi.push_back(Phifrac);
        }
        // NOTE: Calculate Chi-Square
        std::vector<int> ProfileCounts;
        ProfileCounts = MyFun::Hist1D(Phi, bin_cs);
        double chi2_i;
        double Ei = TimeVector.size()/(double) bin_cs;
        for (int i=0; i<ProfileCounts.size(); i++)
        {
            ProfileCounts[i] = std::pow((ProfileCounts[i] - Ei),2)/Ei;
            //Calculate chisquare test, profile counts is not actual profile counts here
        }
        for (auto& n : ProfileCounts) chi2_i += n;
        Chisquare.push_back(chi2_i);

        // clear Vector
        Phi.clear();
    }
    std::cout << std::endl << "FINISH Frequency Search" << std::endl;


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

void CalChisquare() {}

/*void Histogram1D(double* data, double*& double* binsize)
{
}
*/






