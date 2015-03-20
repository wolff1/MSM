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

	Simulator->NumSimulations = Simulator->NumDomains * Simulator->NumMethods;
	Simulator->Methods = (METHOD**) dynmem(Simulator->NumSimulations*sizeof(METHOD*));
	Simulator->Simulations = (SIMULATION**) dynmem(Simulator->NumSimulations*sizeof(SIMULATION*));
	SimulationMatrix = (double**) dynarr_d(Simulator->NumDomains, Simulator->NumMethods);
	NumMethods = 0;
	NumSims = 0;
	for (i = 0; i < Simulator->NumMethods; i++)
	{
		printf("Which method type (NAIVE=0,MSM=1) is method #%hd: ", i);
		scanf("%hd", &SelectedMethod);
//	FIXME - ADD SOMETHING FOR METHOD OPTIONS
		for (j = 0; j < Simulator->NumDomains; j++)
		{
			printf("Apply method <%hd> to domain <%s>? [0|1]: ", i, Simulator->Domains[j]->Name);
			scanf("%hd", &Answer);
			if (Answer)
			{
				FoundMatch = 0;
				for (k = 0; k < j; k++)
				{
					if (SimulationMatrix[k][i] == Simulator->Domains[j]->Radius)
					{
printf("FOUND A MATCH!\n");
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
	dynfree(SimulationMatrix[0]);
	dynfree(SimulationMatrix);

	Simulator->NumSimulations = NumSims;
	Simulator->NumMethods = NumMethods;
	printf("Simulator->NumSimulations = %hd, Methods = %hd\n", Simulator->NumSimulations, Simulator->NumMethods);

	//	Then, run the simulations
	simulator_run_simulations(Simulator);
}

void simulator_examine_results(SIMULATOR* Simulator)
{
	double		abso = 0.0;
	double		rela = 0.0;
	long		i = 0;
//	printf("These are your results, jerk-face!\n");

	//	at most one naive method per domain (run if no file or file is more than a day old)
	//	for each domain, compare the different methods to naive and each other
	//		energy errors (abs / rel)
	//		force errors (abs / rel)
/*
	Simulator->Simulations[0]->Domain->Particles->U;	//	Exact ENERGY value, scalar
	Simulator->Simulations[0]->Domain->Particles->f;	//	Exact FORCE vector, Nx3 vector
	Simulator->Simulations[1]->Domain->Particles->U;	//	MSM O(n) ENERGY value, scalar
	Simulator->Simulations[1]->Domain->Particles->f;	//	MSM O(n) FORCE vector, Nx3 vector
	Simulator->Simulations[2]->Domain->Particles->U;	//	MSM O(nlogn) ENERGY value, scalar
	Simulator->Simulations[2]->Domain->Particles->f;	//	MSM O(nlogn) FORCE vector, Nx3 vector

	printf("%e, %e, %e\n",	Simulator->Simulations[0]->Domain->Particles->U,
							Simulator->Simulations[1]->Domain->Particles->U,
							Simulator->Simulations[2]->Domain->Particles->U);
*/
	printf("0: %32.16f\n",	Simulator->Simulations[0]->Domain->Particles->U);
	abso = fabs(Simulator->Simulations[1]->Domain->Particles->U - Simulator->Simulations[0]->Domain->Particles->U);
	rela = abso / fabs(Simulator->Simulations[0]->Domain->Particles->U);
	printf("ENERGY: abs O(n): %e, rel O(n): %e\n", abso, rela);
	abso = 0.0;
	for (i = 0; i < Simulator->Simulations[0]->Domain->Particles->N; i++)
	{
		abso = MAX(abso,
					fabs(Simulator->Simulations[0]->Domain->Particles->f[i][0] - Simulator->Simulations[1]->Domain->Particles->f[i][0]));
		abso = MAX(abso,
					fabs(Simulator->Simulations[0]->Domain->Particles->f[i][1] - Simulator->Simulations[1]->Domain->Particles->f[i][1]));
		abso = MAX(abso,
					fabs(Simulator->Simulations[0]->Domain->Particles->f[i][2] - Simulator->Simulations[1]->Domain->Particles->f[i][2]));
	}
	printf("FORCE: MAX abs O(n): %e, rel O(n): %e\n", abso, -1.0);
	printf("\n");

	printf("1: %32.16f\n",	Simulator->Simulations[1]->Domain->Particles->U);
	abso = fabs(Simulator->Simulations[2]->Domain->Particles->U - Simulator->Simulations[0]->Domain->Particles->U);
	rela = abso / fabs(Simulator->Simulations[0]->Domain->Particles->U);
	printf("ENERGY: abs O(nlogn): %e, rel O(n): %e\n", abso, rela);
	abso = 0.0;
	for (i = 0; i < Simulator->Simulations[0]->Domain->Particles->N; i++)
	{
		abso = MAX(abso,
					fabs(Simulator->Simulations[0]->Domain->Particles->f[i][0] - Simulator->Simulations[2]->Domain->Particles->f[i][0]));
		abso = MAX(abso,
					fabs(Simulator->Simulations[0]->Domain->Particles->f[i][1] - Simulator->Simulations[2]->Domain->Particles->f[i][1]));
		abso = MAX(abso,
					fabs(Simulator->Simulations[0]->Domain->Particles->f[i][2] - Simulator->Simulations[2]->Domain->Particles->f[i][2]));
	}
	printf("FORCE: MAX abs O(nlong): %e, rel O(n): %e\n", abso, -1.0);
	printf("\n");

	printf("2: %32.16f\n",	Simulator->Simulations[2]->Domain->Particles->U);
	abso = fabs(Simulator->Simulations[1]->Domain->Particles->U - Simulator->Simulations[2]->Domain->Particles->U);
	rela = abso / fabs(Simulator->Simulations[0]->Domain->Particles->U);
	printf("ENERGY: abs O(n) vs O(nlogn): %e, rel: %e\n", abso, rela);
	abso = 0.0;
	for (i = 0; i < Simulator->Simulations[0]->Domain->Particles->N; i++)
	{
		abso = MAX(abso,
					fabs(Simulator->Simulations[1]->Domain->Particles->f[i][0] - Simulator->Simulations[2]->Domain->Particles->f[i][0]));
		abso = MAX(abso,
					fabs(Simulator->Simulations[1]->Domain->Particles->f[i][1] - Simulator->Simulations[2]->Domain->Particles->f[i][1]));
		abso = MAX(abso,
					fabs(Simulator->Simulations[1]->Domain->Particles->f[i][2] - Simulator->Simulations[2]->Domain->Particles->f[i][2]));
	}
	printf("FORCE: MAX abs O(n) vs O(nlogn): %e, rel: %e\n", abso, -1.0);
	printf("\n");
}

void simulator_uninitialize(SIMULATOR* Simulator)
{
	//	Free dynamically allocated memorys
	short		i = 0;

	for (i = 0; i < Simulator->NumDomains; i++)
	{
//printf("Simulator uninitializing Domain <%hd>\n", i);
		simulation_domain_uninitialize(Simulator->Domains[i]);
		dynfree(Simulator->Domains[i]);
	}
//printf("Simulator uninitializing Domains list\n");
	dynfree(Simulator->Domains);
	Simulator->NumDomains = 0;

	for (i = 0; i < Simulator->NumMethods; i++)
	{
//printf("Simulator uninitializing Method <%hd>\n", i);
		method_uninitialize(Simulator->Methods[i]);
		dynfree(Simulator->Methods[i]);
	}
//printf("Simulator uninitializing Methods list\n");
	dynfree(Simulator->Methods);
	Simulator->NumMethods = 0;

	for (i = 0; i < Simulator->NumSimulations; i++)
	{
//printf("Simulator uninitializing Simulation <%hd>\n", i);
		simulation_uninitialize(Simulator->Simulations[i]);
		dynfree(Simulator->Simulations[i]);
	}
//printf("Simulator uninitializing Simulations list\n");
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
	printf("\n");
}

void simulator_run_simulation(SIMULATOR* Simulator, short Index)
{
	//	Start a simulation in its own thread
	simulation_run(Simulator->Simulations[Index]);
}

//	End of file