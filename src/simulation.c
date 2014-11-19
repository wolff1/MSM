//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulation.c - 
*/

#include "simulation.h"

//	EXTERNAL Methods
void simulation_initialize(SIMULATION* Sim)
{
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain
}

void simulation_run(SIMULATION* Sim)
{
	//	Compute forces and electrostatic energy of domain using method
}

void simulation_uninitialize(SIMULATION* Sim)
{
	//	Free dynamically allocated memory
}

//	INTERNAL Methods
void simulation_step(SIMULATION* Sim)
{
	//	Use the calculated forces / momentum information to move the particles
}

//	End of file