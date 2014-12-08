//-------|---------|---------|---------|---------|---------|---------|---------|
/*
particle_collection.h - 
*/

#ifndef PARTICLE_COLLECTION_H
#define PARTICLE_COLLECTION_H

#include "all.h"
#include "memory.h"
#include "particle.h"

typedef struct
{
	long			N;				//	Number of particles
	double			U;				//	Energy induced by particles
	PARTICLE*		r;				//	Positions and charges of particles
	double**		f;				//	Per particle induced force
	double*			m;				//	Masses of particles
	double**		v;				//	Per particle velocity
	double			UnitConverter;	//	Charge scalar for unit conversion
} PARTICLE_COLLECTION;

//	EXTERNAL Methods
void particle_collection_initialize(PARTICLE_COLLECTION* Pc, long N, double UnitConverter);
void particle_collection_move(PARTICLE_COLLECTION* Pc);
void particle_collection_copy(PARTICLE_COLLECTION* DstParticles, PARTICLE_COLLECTION* SrcParticles);
void particle_collection_uninitialize(PARTICLE_COLLECTION* Pc);

#endif

//	End of file