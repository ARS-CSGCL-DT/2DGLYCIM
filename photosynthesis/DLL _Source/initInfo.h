#pragma once
#ifndef _INITINFO_H_
#define _INITINFO_H_
#define MINUTESPERDAY (24*60);

// Contains field and  cultivar information

struct TInitInfo
{
public:
	TInitInfo()
	{

	}

	double latitude, longitude, altitude;// longitude need hard coded
	
	double NRATIO; 
	
};

#endif