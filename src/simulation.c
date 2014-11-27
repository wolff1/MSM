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
	METHOD*					TmpMethod = NULL;

	assert(Simulation != NULL);
	assert(Domain != NULL);
	assert(Method != NULL);

	Simulation->Id = Id;
	Simulation->TimeSteps = TimeSteps;

	//	Make *separate* copies of Domain and Method
	Simulation->Domain = (SIMULATION_DOMAIN*) dynmem(sizeof(SIMULATION_DOMAIN));
	simulation_domain_copy(Domain, Simulation->Domain);

/*
	TmpMethod = (METHOD*) dynmem(sizeof(METHOD));
	method_copy(Method, TmpMethod);
	(*Simulation)->Method = TmpMethod;
*/
}

void simulation_run(SIMULATION* Simulation)
{
	//	Compute forces and electrostatic energy of domain using method
	long		i = 0;

	for (i = 0; i < Simulation->TimeSteps; i++)
	{
//		Simulation->Method->evaluate(Simulation->Method);
		simulation_step(Simulation);

		if (0/*Domain enlarged*/)
		{
			//	Resize K
			Simulation->Method->preprocess(Simulation->Method, Simulation->Domain->Radius);
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

//	method_uninitialize(Simulation->Method);
}

//	INTERNAL Methods
void simulation_step(SIMULATION* Simulation)
{
	//	Use the calculated forces / momentum information to move the particles
}

//	End of file
