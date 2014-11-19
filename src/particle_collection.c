//-------|---------|---------|---------|---------|---------|---------|---------|
/*
particle_collection.c - 
*/

#include "particle_collection.h"

//	EXTERNAL Methods
void particle_collection_initialize(PARTICLE_COLLECTION** Pc)
{
	//	A particle collection is a list of particles and any relevant properties (mass, velocity, etc)
}

void particle_collection_move(PARTICLE_COLLECTION* Pc)
{
	//	This method will update the particle positions at each time step
}

void particle_collection_uninitialize(PARTICLE_COLLECTION* Pc)
{
	//	Free dynamically allocated memory
}

//	End of file