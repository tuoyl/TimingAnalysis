#include <iostream>
#include <string.h>
#include "fitsio.h"
#include <vector>
#include <cmath>
#include "/Users/tuoyouli/cpplibrary/MyFun.h"

//prototype
void PrintError(int status);
//

/* global variables */
/* TODO: read parameters from parfile */ 
double MinFreq = 26;
double MaxFreq = 28;
double FreqStep= 1e-7;
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
    for (int i=0; i<FreqArray.size(); i++)
    {
        for (int i=0; i<TimeVector.size(); i++)
        {
            double Phitmp = (TimeVector[i] - PEpoch)*FreqArray[i] +
                    0.5*std::pow((TimeVector[i] - PEpoch), 2)*F1 +
                    (1.0/6.0)*std::pow((TimeVector[i] - PEpoch), 3)*F2;
            double Phifrac, Phiint;
            Phifrac = std::modf(Phitmp, &Phiint);
            std::cout << "frac" << Phifrac << std::endl;
            Phi.push_back(Phifrac);
        }
        std::vector<int> p_num;
        p_num = MyFun::Hist1D(Phi, 20);
        double chi2tmp;
        for (auto& n : p_num) chi2tmp += n;
        Chisquare.push_back(chi2tmp);
        printf("%f \n", Chisquare);

        // clear Vector
        Phi.clear();
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

/*void Histogram1D(double* data, double*& double* binsize)
{
}
*/






