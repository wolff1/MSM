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

typedef struct
{
	SIMULATION_DOMAIN*		Domain;
	METHOD*					Method;
	SIMULATION*				Simulation;
	long					TimeSteps;
//	OUTPUT*					Output;
} SIMULATOR;

//	EXTERNAL Methods
void simulator_initialize(SIMULATOR* Sim);
void simulator_run(SIMULATOR* Sim);
void simulator_run_all(SIMULATOR* Sim);
void simulator_uninitialize(SIMULATOR* Sim);

#endif

//	End of file