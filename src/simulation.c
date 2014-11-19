//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulation.c - 
*/

#include "simulation.h"

//	EXTERNAL Methods
void simulation_initialize(SIMULATION** Simulation, SIMULATION_DOMAIN* Domain, METHOD* Method)
{
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain
	assert(*Simulation == NULL);
	assert(Domain != NULL);
	assert(Method != NULL);
	*Simulation = (SIMULATION*) dynmem(sizeof(SIMULATION));
	(*Simulation)->TimeSteps = 1;
	Method->preprocess((void*)Method);
}

void simulation_run(SIMULATION* Simulation)
{
	//	Compute forces and electrostatic energy of domain using method
	long		i = 0;

	for (i = 0; i < Simulation->TimeSteps; i++)
	{
		if (0/* OR Domain enlarged*/)
		{
			Simulation->Method.preprocess(&Simulation->Method);
		}

		//Simulation->Method.evaluate(&Simulation->Method);
		simulation_step(Simulation);
	}
	//Sim->Method.uninitialize(&Sim->Method); //<--- happens in simulator
}

void simulation_uninitialize(SIMULATION* Simulation)
{
	//	Free dynamically allocated memory
	dynfree(Simulation);
}

//	INTERNAL Methods
void simulation_step(SIMULATION* Simulation)
{
	//	Use the calculated forces / momentum information to move the particles
}

//	End of file