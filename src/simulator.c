//-------|---------|---------|---------|---------|---------|---------|---------|
/*
simulator.c - 
*/

#include "simulator.h"

//	EXTERNAL Methods
void simulator_initialize(SIMULATOR* Simulator)
{
	//	A simulator object contains a) simulations, b) domains, c) methods, and d) other stuff
	//	A simulation is the repeated application of a specific method to a certain domain
	//	A domain is a system of particles, and any other relevant information
	//	A method computes the forces and electrostatic energy within a domain
	assert(Simulator != NULL);

	Simulator->NumDomains = 0;
	Simulator->Domains = NULL;

	Simulator->NumMethods = 0;
	Simulator->Methods = NULL;

	Simulator->NumSimulations = 0;
	Simulator->Simulations = NULL;
}

void simulator_run(SIMULATOR* Simulator)
{
	//	Build list of simulations to run
	short					i = 0;
	short					j = 0;
	size_t					MethodSize = 0;
	void*					Init = NULL;
	SIMULATION_DOMAIN*		SimulationDomain = NULL;
	METHOD*					Method = NULL;
	SIMULATION*				Simulation = NULL;
	short					DomainIdx = 0;
	short					MethodIdx = 0;
	char					DomainFileName[128] = {0};
	short					SelectedMethod = 0;

	short					NumSims = 0;
	short					NumMethods = 0;
	double**				SimulationMatrix = NULL;
	short					Answer = 0;
	short					k = 0;
	short					FoundMatch = 0;

	printf("How many domains? ");
	scanf("%hd", &Simulator->NumDomains);
	Simulator->Domains = (SIMULATION_DOMAIN**) dynmem(Simulator->NumDomains*sizeof(SIMULATION_DOMAIN*));
	for (i = 0; i < Simulator->NumDomains; i++)
	{
		//	Which domain? (consider filename to be unique)
		//		Create and initialize domain?
		printf("Enter the file name for domain #%hd: ", i);
		scanf("%s", DomainFileName);

		SimulationDomain = (SIMULATION_DOMAIN*) dynmem(sizeof(SIMULATION_DOMAIN));
		simulation_domain_initialize(SimulationDomain, i, DomainFileName);
		Simulator->Domains[i] = SimulationDomain;
	}

	printf("How many methods? ");
	scanf("%hd", &Simulator->NumMethods);
	Simulator->Methods = (METHOD**) dynmem(Simulator->NumMethods*Simulator->NumDomains*sizeof(METHOD*));
	NumMethods = 0;
/*
	for (i = 0; i < Simulator->NumMethods; i++)
	{
		printf("Please enter which method (NAIVE=0,MSM=1) for method #%hd: ", i);
		scanf("%hd", &SelectedMethod);

		//	Which method? (consider combination of parameters to be unique)
		//		Create and initialize method?
		switch (SelectedMethod)
		{
		case 1:		//	MSM
			Init = &msm_initialize;
			MethodSize = sizeof(MSM);
			break;
		default:	//	NAIVE
			Init = &naive_initialize;
			MethodSize = sizeof(NAIVE);
		}
		Method = (METHOD*) dynmem(MethodSize);
		method_initialize((void*)Method, Init, i);
		Simulator->Methods[i] = Method;
		NumMethods++;
	}
*/
	//	Assume each method will be used for each domain
	Simulator->NumSimulations = Simulator->NumDomains * Simulator->NumMethods;
	Simulator->Simulations = (SIMULATION**) dynmem(Simulator->NumSimulations*sizeof(SIMULATION*));
	SimulationMatrix = (double**) dynarr_d(Simulator->NumDomains, Simulator->NumMethods);
	NumSims = 0;
	for (i = 0; i < Simulator->NumMethods; i++)
	{
		printf("Please enter which method (NAIVE=0,MSM=1) for method #%hd: ", i);
		scanf("%hd", &SelectedMethod);

		for (j = 0; j < Simulator->NumDomains; j++)
		{
			//	Apply method i to domain j? If yes:
			//		if method i NOT initialized for domain[j].size
			//			method_preprocess
			//		simulation = method + domain
			//		initialize simulation
			printf("Apply method <%hd> to domain <%s>? [0|1]: ", i, Simulator->Domains[j]->Name);
			scanf("%hd", &Answer);
			if (Answer)
			{
				FoundMatch = 0;
				for (k = 0; k < j; k++)
				{
					if (SimulationMatrix[k][i] == Simulator->Domains[j]->Radius)
					{
						//	We found another domain with the same size for the same method!
						FoundMatch = 1;
						break;
					}
				}

				if (!FoundMatch)
				{
					//	Allocate a new Method because the one we need does not exist
					switch (SelectedMethod)
					{
					case 1:		//	MSM
						Init = &msm_initialize;
						MethodSize = sizeof(MSM);
						break;
					default:	//	NAIVE
						Init = &naive_initialize;
						MethodSize = sizeof(NAIVE);
					}
					Method = (METHOD*) dynmem(MethodSize);
					method_initialize((void*)Method, Init, NumMethods);
					Simulator->Methods[NumMethods] = Method;

					//	None of the previous domains have the same size as the current one
					Simulator->Methods[NumMethods]->preprocess(Simulator->Methods[NumMethods], Simulator->Domains[j]->Radius);
					SimulationMatrix[j][i] = Simulator->Domains[j]->Radius;

					NumMethods++;
				}
				else
				{
					//	Use the pre-existing method for this simulation
					//		-> Simulator->Methods[i] has already been pre-processed for the correct size
				}

				//	Create simulation (these could happen in parallel with OpenMP)
				Simulation = (SIMULATION*) dynmem(sizeof(SIMULATION));
				simulation_initialize(Simulation, Simulator->Domains[j], Simulator->Methods[NumMethods-1], NumSims, 1);	//FIXME <-- 1 is TimeSteps
				Simulator->Simulations[NumSims] = Simulation;
				NumSims++;
			}
		}
	}
	Simulator->NumSimulations = NumSims;
	Simulator->NumMethods = NumMethods;
printf("Simulator->NumSimulations = %ld, Methods = %ld\n", Simulator->NumSimulations, Simulator->NumMethods);
	dynfree(SimulationMatrix[0]);
	dynfree(SimulationMatrix);
/*
	//	Assume each method will be used for each domain
	Simulator->NumSimulations = Simulator->NumDomains * Simulator->NumMethods;
	Simulator->Simulations = (SIMULATION**) dynmem(Simulator->NumSimulations*sizeof(SIMULATION*));
	for (i = 0; i < Simulator->NumSimulations; i++)
	{
		DomainIdx = i / Simulator->NumMethods;
		MethodIdx = i % Simulator->NumMethods;

		//	Create simulation (these could happen in parallel with OpenMP)
		Simulation = (SIMULATION*) dynmem(sizeof(SIMULATION));
		simulation_initialize(Simulation, Simulator->Domains[DomainIdx], Simulator->Methods[MethodIdx], i, 1);	//FIXME <-- 1 is TimeSteps
		Simulator->Simulations[i] = Simulation;
	}
*/
	//	Then, run the simulations
	simulator_run_simulations(Simulator);
}

void simulator_uninitialize(SIMULATOR* Simulator)
{
	//	Free dynamically allocated memorys
	short		i = 0;

	for (i = 0; i < Simulator->NumDomains; i++)
	{
printf("Simulator uninitializing Domain <%hd>\n", i);
		simulation_domain_uninitialize(Simulator->Domains[i]);
		dynfree(Simulator->Domains[i]);
	}
printf("Simulator uninitializing Domains list\n");
	dynfree(Simulator->Domains);
	Simulator->NumDomains = 0;

	for (i = 0; i < Simulator->NumMethods; i++)
	{
printf("Simulator uninitializing Method <%hd>\n", i);
		method_uninitialize(Simulator->Methods[i]);
		dynfree(Simulator->Methods[i]);
	}
printf("Simulator uninitializing Methods list\n");
	dynfree(Simulator->Methods);
	Simulator->NumMethods = 0;

	for (i = 0; i < Simulator->NumSimulations; i++)
	{
printf("Simulator uninitializing Simulation <%hd>\n", i);
		simulation_uninitialize(Simulator->Simulations[i]);
		dynfree(Simulator->Simulations[i]);
	}
printf("Simulator uninitializing Simulations list\n");
	dynfree(Simulator->Simulations);
	Simulator->NumSimulations = 0;
}

//	INTERNAL Methods
void simulator_run_simulations(SIMULATOR* Simulator)
{
	//	Start all simulations within the simulator
	//		-> Use OpenMP to multi-thread simulations at this level?
	short		i = 0;

	for (i = 0; i < Simulator->NumSimulations; i++)
	{
printf("Simulator running Simulation <%hd>\n", i);
		simulator_run_simulation(Simulator, i);
	}
}

void simulator_run_simulation(SIMULATOR* Simulator, short Index)
{
	//	Start a simulation in its own thread
	simulation_run(Simulator->Simulations[Index]);
}

//	End of file