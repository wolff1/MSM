//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulation.c - 
*/

#include "simulation.h"

//	EXTERNAL Methods
void simulation_initialize(SIMULATION* Simulation, SIMULATION_DOMAIN* Domain, METHOD* Method, short Id, long TimeSteps)
{
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain

	assert(Simulation != NULL);
	assert(Domain != NULL);
	assert(Method != NULL);

	Simulation->Id = Id;
	Simulation->TimeSteps = TimeSteps;

	//	Make *separate* copy of Domain
	Simulation->Domain = (SIMULATION_DOMAIN*) dynmem(sizeof(SIMULATION_DOMAIN));
	simulation_domain_copy(Simulation->Domain, Domain);

	//	Make *separate* copy of Method
	Simulation->Method = (METHOD*) dynmem(Method->Size);
	method_copy(Simulation->Method, Method);
}

void simulation_run(SIMULATION* Simulation)
{
	//	Compute forces and electrostatic energy of domain using method
	long		i = 0;

	for (i = 0; i < Simulation->TimeSteps; i++)
	{
		//	Have <METHOD> evaluate the energy and forces for the <DOMAIN>
		(*Simulation->Method->evaluate)(Simulation->Method, Simulation->Domain->Particles->N, Simulation->Domain->Particles->r);

		//	Handle time integration
		simulation_step(Simulation);

		if (0/*Domain enlarged*/)
		{
			//	FIXME - Probably should have SIMULATION_DOMAIN do the resizing
			//	Resize K
			(*Simulation->Method->preprocess)(Simulation->Method, Simulation->Domain->Radius);
		}
	}
}

void simulation_uninitialize(SIMULATION* Simulation)
{
	//	Free dynamically allocated memory
printf("before simulation_domain_uninit\n");
  	simulation_domain_uninitialize(Simulation->Domain);
printf("before dynfree(Simulation->Domain)\n");
	dynfree(Simulation->Domain);

printf("before method_uninit\n");
	method_uninitialize(Simulation->Method);
printf("before dynfree(Simulation->Method)\n");
	dynfree(Simulation->Method);
}

//	INTERNAL Methods
void simulation_step(SIMULATION* Simulation)
{
	//	Have SIMULATION_DOMAIN Use the calculated forces / momentum information to move the particles
}

//	End of file
