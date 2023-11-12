#pragma once

#ifndef _CONTROLLER_H_
#define _CONTROLLER_H_

#include "weather.h"
#include "initinfo.h"
#include "gas exchange new.h"


class CController
{

private:

	TInitInfo			initInfo;

	TWeather* weather;

	CGas_exchange_new* sunlit;
	CGas_exchange_new* shaded;

	double sunLAI, shadeLAI;
	double photosynthesis_gross; // gross photosynthesis, umolCO2 m-2 ground s-1
	double photosynthesis_grosssunlitleaf;//pg umol co2 m-2 leaf s-1, sunlit leaves
	double photosynthesis_grossshadedleaf;//pg umol co2 m-2 leaf s-1, shaded leaves
	double photosynthesis_net; // gross photosynthesis, umolCO2 m-2 ground s-1	
	double transpiration; //instanteous transpiration, mmol H2O m-2 ground s-1
	double sunPg, shadePg; //store Pg (umol co2 m-2 leaf s-1) estimated for all sunlit and shaded leaves
	double sunPFD, shadePFD;//instantaneous incident, shaded and sunlit PAR values in umol PAR m-2 s-1, SRAD in W m-2
	double internal_CO2; //hourly leaf sub-stomatal CO2 concetration (ppm)
	double photosynthesis_netsunlitleaf; //net p.s., umol co2 m-2 leaf s-1
	double photosynthesis_netshadedleaf; //net p.s., umol co2 m-2 leaf s-1
	double temperature; // leaf T, hourly
	double transpiration_sunlitleaf; //leaf et, mmol H2O m-2 leaf s-1
	double transpiration_shadedleaf; //leaf et, mmol h2o m-2 leaf s-1
	double conductance; //weighted stomatal conductance from gas exchange class for current hour, mmol H2o m-2 s-1
	double VPD;
	double psistress_gs_factor;
	double sunlitLAI, shadedLAI; // sunlit and shaded LAI values
	double sunlitPFD, shadedPFD;
	double TLAI;
	double temp1;
	double Ags;
public:

	CController();
	CController(const TInitInfo);
	~CController();
	

	int run(const TWeather &, TInitInfo);
	CGas_exchange_new* get_sunlit() { return this->sunlit; }
	CGas_exchange_new* get_shaded() { return this->shaded; } // get access to the pointers that point to the sunlit/shade leaves 

	void set_sunPg(double x) { sunPg = x; }
	void set_shadePg(double x) { shadePg = x; }
	double get_photosynthesis_gross() { return photosynthesis_gross; }
	double get_photosynthesis_net() { return photosynthesis_net; }
	double get_transpiration() { return transpiration; }
	double get_temperature() { return temperature; } // leaf temperature
	double get_TLAI() { return TLAI; }
	double get_photosynthesis_netsunlitleaf() { return photosynthesis_netsunlitleaf; }
	double get_photosynthesis_netshadedleaf() { return photosynthesis_netshadedleaf; }
	double get_transpiration_sunlitleaf() { return transpiration_sunlitleaf; }
	double get_transpiration_shadedleaf() { return transpiration_shadedleaf; }
	double get_sunlitLAI() { return sunlitLAI; }
	double get_shadedLAI() { return shadedLAI; }
	double get_temp1() { return temp1; }
	double get_Ags() { return Ags; }
	};
#endif

