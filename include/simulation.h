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
	short					Id;
	SIMULATION_DOMAIN*		Domain;
	METHOD*					Method;
	long					TimeSteps;
//	OUTPUT					Output;
} SIMULATION;

//	EXTERNAL Methods
void simulation_initialize(SIMULATION** Simulation, SIMULATION_DOMAIN* Domain, METHOD* Method, short Id, long TimeSteps);
void simulation_run(SIMULATION* Simulation);
void simulation_uninitialize(SIMULATION* Simulation);

//	INTERNAL Methods
void simulation_step(SIMULATION* Simulation);

#endif

//	End of file