//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulation.c - 
*/

#include "simulation.h"

//	EXTERNAL Methods
void simulation_initialize(SIMULATION** Simulation, SIMULATION_DOMAIN* Domain, METHOD* Method, short Id, long TimeSteps)
{
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain
	assert(*Simulation == NULL);
	assert(Domain != NULL);
	assert(Method != NULL);

	(*Simulation) = (SIMULATION*) dynmem(sizeof(SIMULATION));

	(*Simulation)->Id = Id;
	(*Simulation)->TimeSteps = TimeSteps;

	//	Make *separate* copies of Domain and Method
	//		-> They've already been initialized (i.e. preprocessing, etc)
//	domain_copy(Domain, &Simulation->Domain);
//	method_copy(Method, &Simulation->Method);
}

void simulation_run(SIMULATION* Simulation)
{
	//	Compute forces and electrostatic energy of domain using method
	long		i = 0;

	for (i = 0; i < Simulation->TimeSteps; i++)
	{
		if (0/*Domain enlarged*/)
		{
			Simulation->Method->preprocess(Simulation->Method);
		}

//		Simulation->Method->evaluate(Simulation->Method);
		simulation_step(Simulation);
	}
}

void simulation_uninitialize(SIMULATION* Simulation)
{
	//	Free dynamically allocated memory
//	simulation_domain_uninitialize(Simulation->Domain);
//	method_uninitialize(Simulation->Method);
	dynfree(Simulation);
}

//	INTERNAL Methods
void simulation_step(SIMULATION* Simulation)
{
	//	Use the calculated forces / momentum information to move the particles
}

//	End of file