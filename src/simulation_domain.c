//-------|---------|---------|---------|---------|---------|---------|---------|
/*
domain.c - 
*/

#include "simulation_domain.h"

//	EXTERNAL Methods
void simulation_domain_initialize(SIMULATION_DOMAIN* Domain);
void simulation_domain_compute_energy(SIMULATION_DOMAIN* Domain);
void simulation_domain_compute_forces(SIMULATION_DOMAIN* Domain);
void simulation_domain_uninitialize(SIMULATION_DOMAIN* Domain);

//	End of file