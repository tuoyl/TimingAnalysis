#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

void ReadParFile(std::string ParfileName, double& PEPOCH, double& F0, double& FreqRange, double& FreqStep, int& BinNum)
{
    std::ifstream inFile;
    std::string STRING;
    std::string StringParValue;
    double ParValue;
    inFile.open(ParfileName);
    if (!inFile)
    {
        std::cout << "Unable to open file " << ParfileName << std::endl;
        exit(1);
    }

    std::cout << "########## Parameters ############" << std::endl;
    while (inFile >> STRING)
    {
        if (STRING == "NAME")
        {
            inFile >> StringParValue;
            std::cout << STRING << ":" << StringParValue << std::endl;
        }
        if (STRING == "RAJ")
        {
            inFile >> StringParValue;
            std::cout << STRING << ":" << StringParValue << std::endl;
        }
        if (STRING == "DECJ")
        {
            inFile >> StringParValue;
            std::cout << STRING << ":" << StringParValue << std::endl;
        }
        if (STRING == "PEPOCH")
        {
            inFile >> ParValue;
            std::cout << STRING << ":" << ParValue << std::endl;
            PEPOCH = ParValue;
        }
        if (STRING == "F0")
        {   
            inFile >> ParValue;
            std::cout << STRING << ":" << ParValue << std::endl;
            F0 = ParValue;
        }
        if (STRING == "SEARCH_RANGE")
        {   
            inFile >> ParValue;
            std::cout << STRING << ":" << ParValue << std::endl;
            FreqRange = ParValue;
        }
        if (STRING == "SEARCH_STEP")
        {   
            inFile >> ParValue;
            std::cout << STRING << ":" << ParValue << std::endl;
            FreqStep = ParValue;
        }
        if (STRING == "BinNumber")
        {   
            inFile >> ParValue;
            std::cout << STRING << ":" << ParValue << std::endl;
            BinNum = ParValue;
        }
    }
}
