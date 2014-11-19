//-------|---------|---------|---------|---------|---------|---------|---------|
/*
particle_collection.c - 
*/

#include "particle_collection.h"

//	EXTERNAL Methods
void particle_collection_initialize(PARTICLE_COLLECTION** Pc, long N, double UnitConverter)
{
	//	A particle collection is a list of particles and any relevant properties (mass, velocity, etc)
	assert(*Pc == NULL);

	*Pc = (PARTICLE_COLLECTION*) dynmem(sizeof(PARTICLE_COLLECTION));

	(*Pc)->N = N;
	(*Pc)->r = (PARTICLE*) dynmem((*Pc)->N*sizeof(PARTICLE));
	(*Pc)->m = (double*) dynvec((*Pc)->N, sizeof(double));
	(*Pc)->v = (double**) dynarr_d((*Pc)->N, 3);
	(*Pc)->UnitConverter = UnitConverter;
}

void particle_collection_move(PARTICLE_COLLECTION* Pc)
{
	//	This method will update the particle positions at each time step
}

void particle_collection_uninitialize(PARTICLE_COLLECTION* Pc)
{
	//	Free dynamically allocated memory
	dynfree(Pc->v[0]);
	dynfree(Pc->v);
	dynfree(Pc->m);
	dynfree(Pc->r);
	dynfree(Pc);
}

//	End of file