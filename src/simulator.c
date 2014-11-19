//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulator.c - 
*/

#include "simulator.h"

//	EXTERNAL Methods
void simulator_initialize(SIMULATOR** Simulator)
{
	//	A simulator object contains a) simulations, b) domains, c) methods, and d) other stuff
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain
	assert(*Simulator == NULL);

	(*Simulator) = (SIMULATOR*) dynmem(sizeof(SIMULATOR));

	(*Simulator)->NumDomains = 0;
	(*Simulator)->Domains = NULL;

	(*Simulator)->NumMethods = 0;
	(*Simulator)->Methods = NULL;

	(*Simulator)->NumSimulations = 0;
	(*Simulator)->Simulations = NULL;
}

void simulator_run(SIMULATOR* Simulator)
{
	//	Build list of simulations to run
	short					i = 0;
	size_t					MethodSize = 0;
	void*					Init = NULL;
	SIMULATION_DOMAIN*		SimulationDomain = NULL;
	METHOD*					Method = NULL;
	SIMULATION*				Simulation = NULL;
	short					DomainIdx = 0;
	short					MethodIdx = 0;
	char					DomainFileName[128] = {0};
	short					SelectedMethod = 0;

	printf("How many domains? ");
	scanf("%hd", &Simulator->NumDomains);
	Simulator->Domains = (SIMULATION_DOMAIN**) dynmem(Simulator->NumDomains*sizeof(SIMULATION_DOMAIN*));
	for (i = 0; i < Simulator->NumDomains; i++)
	{
		//	Which domain? (consider filename to be unique)
		//		Create and initialize domain?
		printf("Enter the file name for domain #%hd: ", i);
		scanf("%s", DomainFileName);

		SimulationDomain = NULL;
		simulation_domain_initialize(&SimulationDomain, i, DomainFileName);
		Simulator->Domains[i] = SimulationDomain;
	}

	printf("How many methods? ");
	scanf("%hd", &Simulator->NumMethods);
	Simulator->Methods = (METHOD**) dynmem(Simulator->NumMethods*sizeof(METHOD*));
	for (i = 0; i < Simulator->NumMethods; i++)
	{
		printf("Please enter which method (NAIVE=0,MSM=1) for method #%hd: ", i);
		scanf("%hd", &SelectedMethod);

		//	Which method? (consider combination of parameters to be unique)
		//		Create and initialize method?
		if (SelectedMethod == 1)	//	MSM
		{
			Init = &msm_initialize;
			MethodSize = sizeof(MSM);
		}
		else						//	NAIVE
		{
			Init = &naive_initialize;
			MethodSize = sizeof(NAIVE);
		}

		Method = NULL;
		method_initialize((void*)&Method, MethodSize, Init);
		Method->Id = i;
		Method->preprocess(Method);
		Simulator->Methods[i] = Method;
	}

	//	Assume each method will be used for each domain
	Simulator->NumSimulations = Simulator->NumDomains * Simulator->NumMethods;
	Simulator->Simulations = (SIMULATION**) dynmem(Simulator->NumSimulations*sizeof(SIMULATION*));
	for (i = 0; i < Simulator->NumSimulations; i++)
	{
		DomainIdx = i / Simulator->NumMethods;
		MethodIdx = i % Simulator->NumMethods;

		//	Create simulation (these could happen in parallel with OpenMP)
		Simulation = NULL;
		simulation_initialize(&Simulation, Simulator->Domains[DomainIdx], Simulator->Methods[MethodIdx], i, 1);	//FIXME <-- 1 is TimeSteps
		Simulator->Simulations[i] = Simulation;
	}

	//	Then, run the simulations
	simulator_run_simulations(Simulator);
}

void simulator_uninitialize(SIMULATOR* Simulator)
{
	//	Free dynamically allocated memorys
	short		i = 0;

	for (i = 0; i < Simulator->NumDomains; i++)
	{
		simulation_domain_uninitialize(Simulator->Domains[i]);
	}
	dynfree(Simulator->Domains);
	Simulator->NumDomains = 0;

	for (i = 0; i < Simulator->NumMethods; i++)
	{
		method_uninitialize(Simulator->Methods[i]);

	}
	dynfree(Simulator->Methods);
	Simulator->NumMethods = 0;

	for (i = 0; i < Simulator->NumSimulations; i++)
	{
		simulation_uninitialize(Simulator->Simulations[i]);
	}
	dynfree(Simulator->Simulations);
	Simulator->NumSimulations = 0;

	dynfree(Simulator);
}

//	INTERNAL Methods
void simulator_run_simulations(SIMULATOR* Simulator)
{
	//	Start all simulations within the simulator
	//		-> Use OpenMP to multi-thread simulations at this level?
	short		i = 0;

	for (i = 0; i < Simulator->NumSimulations; i++)
	{
		simulator_run_simulation(Simulator, i);
	}
}

void simulator_run_simulation(SIMULATOR* Simulator, short Index)
{
	//	Start a simulation in its own thread
	simulation_run(Simulator->Simulations[Index]);
}

//	End of file