//-------|---------|---------|---------|---------|---------|---------|---------|
/*
domain.c - 
*/

#include "simulation_domain.h"

//	EXTERNAL Methods
void simulation_domain_initialize(SIMULATION_DOMAIN* Domain)
{
	//	A domain is a system of particles, and any other relevant information
}

void simulation_domain_compute_forces(SIMULATION_DOMAIN* Domain)
{
	//	Pass domain info to method and let it calculate the forces and electrostatic energy
}

void simulation_domain_uninitialize(SIMULATION_DOMAIN* Domain)
{
	//	Free dynamically allocated memory
}

//	End of file