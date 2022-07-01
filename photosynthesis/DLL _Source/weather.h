#pragma once
#ifndef _WEATHER_H_
#define _WEATHER_H_

struct TWeather
{
public:
	TWeather()
	{
		DayOfYear = 1, time = 0.0, CO2 = 370.0, airT = 20.0, psil_ = -0.5, PFD = 0.0, solRad = 0.0,
			RH = 50.0, wind = 1.5, LAI = 0.0, ET_supply = 0.0;
	}
	int DayOfYear;
	double time;
	double CO2;
	double airT; // Air temperature
	double psil_; //Leafwater potential, MPa
	double PFD; 
	double solRad; 
	double RH; //Relative humidity in percent
	double wind; //Windspeed m s-1
	double LAI; 
	double ET_supply;
};
#endif

