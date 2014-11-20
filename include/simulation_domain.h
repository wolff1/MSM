//-------|---------|---------|---------|---------|---------|---------|---------|
/*
domain.h - 
*/

#ifndef	SIMULATION_DOMAIN_H
#define	SIMULATION_DOMAIN_H

#include "float.h"

#include "all.h"
#include "memory.h"
#include "particle_collection.h"

#include "file_format.h"

#define	UNITS_COULOMB				332.0636
	// Coulomb constant for electrostatics, units (AMU*(A/fs)^2)*(A/e^2)
#define MAXLEN_DOMAIN_NAME_STR		128

typedef struct
{
	short					Id;
	PARTICLE_COLLECTION*	Particles;
	PARTICLE				MinimumCoordinates;
	PARTICLE				CenterCoordinates;
	PARTICLE				MaximumCoordinates;
	double					Radius;
	char					Name[MAXLEN_DOMAIN_NAME_STR];
} SIMULATION_DOMAIN;

//	EXTERNAL Methods
void simulation_domain_initialize(SIMULATION_DOMAIN* Domain, short Id, char* FileName);
void simulation_domain_compute_forces(SIMULATION_DOMAIN* Domain);
void simulation_domain_uninitialize(SIMULATION_DOMAIN* Domain);

//	INTERNAL Methods:
void simulation_domain_input_particles(SIMULATION_DOMAIN* Domain);

#endif

//	End of file