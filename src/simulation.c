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
		(*Simulation->Method->evaluate)(Simulation->Method, Simulation->Domain);

		//	Handle time integration
		simulation_step(Simulation);

		if (0/*Domain size changed*/)
		{
			//	FIXME - Probably should have SIMULATION_DOMAIN do the resizing
			//	Resize K
			(*Simulation->Method->preprocess)(Simulation->Method, Simulation->Domain->Radius);
		}
	}

	printf("Energy: %+f\n", Simulation->Domain->Particles->U);
	printf("Forces:\n");
	for (i = 0; i < Simulation->Domain->Particles->N; i++)
	{
		printf("%+f\t%+f\t%+f\n", Simulation->Domain->Particles->f[i][0], Simulation->Domain->Particles->f[i][1], Simulation->Domain->Particles->f[i][2]);
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
	simulation_domain_update(Simulation->Domain);
}

//	End of file
