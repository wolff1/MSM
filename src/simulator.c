//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulator.c - 
*/

#include "simulator.h"

//	EXTERNAL Methods
void simulator_initialize(SIMULATOR* Sim)
{
	//	A simulator object contains a) simulations, b) domains, c) methods, and d) other stuff
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain
}

void simulator_run(SIMULATOR* Sim)
{
	//	Start a simulation in its own thread
}

void simulator_run_all(SIMULATOR* Sim)
{
	//	Start all simulations within the simulator
	//		-> Use OpenMP to multi-thread simulations at this level?
}

void simulator_uninitialize(SIMULATOR* Sim)
{
	//	Free dynamically allocated memorys
}

//	End of file