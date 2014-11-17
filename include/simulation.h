//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulation.h - 
*/

#ifndef	SIMULATION_H
#define	SIMULATION_H

#include "all.h"
#include "simulation_domain.h"
#include "method.h"

typedef struct
{
	SIMULATION_DOMAIN		Domain;
	METHOD					Method;
	long					TimeSteps;
//	OUTPUT					Output;
} SIMULATION;

//	EXTERNAL Methods
void simulation_initialize(SIMULATION* Sim);
void simulation_run(SIMULATION* Sim);
void simulation_step(SIMULATION* Sim);
void simulation_uninitialize(SIMULATION* Sim);

#endif

//	End of file