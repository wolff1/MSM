//-------|---------|---------|---------|---------|---------|---------|---------|
/*
domain.h - 
*/

#ifndef	SIMULATION_DOMAIN_H
#define	SIMULATION_DOMAIN_H

#include "all.h"
#include "memory.h"
#include "particle_collection.h"

typedef struct
{
	short					Id;
	PARTICLE_COLLECTION*	Particles;
	PARTICLE				MinimumCoordinates;
	PARTICLE				CenterCoordinates;
	PARTICLE				MaximumCoordinates;
	double					Radius;
} SIMULATION_DOMAIN;

//	EXTERNAL Methods
void simulation_domain_initialize(SIMULATION_DOMAIN** Domain, short Id, char* FileName);
void simulation_domain_compute_forces(SIMULATION_DOMAIN* Domain);
void simulation_domain_uninitialize(SIMULATION_DOMAIN* Domain);

#endif

//	End of file