//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulator.h - 
*/

#ifndef	SIMULATOR_H
#define	SIMULATOR_H

#include "all.h"
#include "simulation_domain.h"
#include "method.h"
#include "simulation.h"

#include "msm.h"
#include "naive.h"

typedef struct
{
	short					NumDomains;
	short					NumMethods;
	short					NumSimulations;
	SIMULATION_DOMAIN**		Domains;
	METHOD**				Methods;
	SIMULATION**			Simulations;
//	OUTPUT*					Output;
} SIMULATOR;

//	EXTERNAL Methods
void simulator_initialize(SIMULATOR* Simulator);
void simulator_run(SIMULATOR* Simulator);
void simulator_uninitialize(SIMULATOR* Simulator);

//	INTERNAL Methods
void simulator_run_simulations(SIMULATOR* Simulator);
void simulator_run_simulation(SIMULATOR* Simulator, short Index);

#endif

//	End of file