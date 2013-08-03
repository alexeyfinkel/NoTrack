//this is a mini script that calculates the number of xtals in each ring.
//it also makes Root explode in a fireball of death, so don't ever use it


#include <stdio.h>
#include <iostream>

int ringCalc()
{
	int ring;
	std::vector<int> nCrystals(8,0);	
	for(int i=30;i<71;i++)
	{
		for(int j=30;j<71;j++)
		{
			ring = (int)sqrt( (i-50)*(i-50)+(j-50)*(j-50) )-11;
			nCrystals[ring]++;
		}
	}
	
	for(int i=0;i<8;i++)
	{
		std::cout<<"Ring <<"<<i<<" has "<<nCrystals[i]<<" xtals."<<std::endl;
	}
	
	return 0;
}
