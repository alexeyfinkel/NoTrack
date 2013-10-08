#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>


int convert()
{
    //int count;
    char line[1000];
    int ix, iy, iz;
    float calConst;//, uncert;
    std::string junk;
    ifstream xmlFile;
    ofstream textFile;
    
    xmlFile.open("Calibration/MC/MC_EE_smearing_2012.txt");
    if(!xmlFile.is_open())
	{
		std::cout<<"Failed to open Smearing consts. file. Existing."<<std::endl;
		return 1;
	}
    textFile.open("Calibration/MC/MC_SmearingConsts.txt");
    if(!xmlFile.is_open())
	{
		std::cout<<"Failed to open the text file. Existing."<<std::endl;
		return 1;
	}
    
    do
    {
        xmlFile.getline(line, 1000);
        sscanf(line,"%*[^\"]\"%d\"%*[^\"]\"%d\"%*[^\"]\"%d\"%*[^\"]\"%f\"",&ix,&iy,&iz,&calConst);
        if( (ix>25)&&(ix<75)&&(iy>25)&&(iy<75))
        {
            textFile<<ix<<"\t"<<iy<<"\t"<<iz<<"\t"<<calConst<<"\t0"<<std::endl;
        }
        //std::cout<<"("<<ix<<","<<iy<<","<<iz<<") "<<calConst<<std::endl;
        //std::cout<<count<<"  "<<line<<std::endl;
    } while(!xmlFile.eof());
    
    return 0;
}