//-------|---------|---------|---------|---------|---------|---------|---------|
/*
particle_collection.h - 
*/

#ifndef PARTICLE_COLLECTION_H
#define PARTICLE_COLLECTION_H

#include "float.h"

#include "all.h"
#include "memory.h"
#include "particle.h"

typedef struct
{
	long			N;				//	Number of particles
	double			U;				//	Energy induced by particles
	double			UnitConverter;	//	Charge scalar for unit conversion
	PARTICLE*		r;				//	Positions and charges of particles
	double*			m;				//	Masses of particles
	double**		f;				//	Per particle induced force
	double**		v;				//	Per particle velocity
} PARTICLE_COLLECTION;

//	EXTERNAL Methods
void particle_collection_initialize(PARTICLE_COLLECTION* Pc, long N, double UnitConverter);
void particle_collection_update(PARTICLE_COLLECTION* Pc, PARTICLE* MinimumCoordinate, PARTICLE* MaximumCoordinate);
void particle_collection_copy(PARTICLE_COLLECTION* DstParticles, PARTICLE_COLLECTION* SrcParticles);
void particle_collection_uninitialize(PARTICLE_COLLECTION* Pc);

#endif

//	End of file