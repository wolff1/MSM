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
	long			N;
	PARTICLE*		r;
	double*			m;
	double**		v;
	double			UnitConverter;
} PARTICLE_COLLECTION;

//	EXTERNAL Methods
void particle_collection_initialize(PARTICLE_COLLECTION** Pc, long N, double UnitConverter);
void particle_collection_move(PARTICLE_COLLECTION* Pc);
void particle_collection_uninitialize(PARTICLE_COLLECTION* Pc);

#endif

//	End of file