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
					method_initialize((void*)Method, Init, NumMethods, NULL, NULL);	//	FIXME: 2 NULLs will not work!
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

void simulator_run_water(SIMULATOR* Simulator)
{
	short					i = 0;
	size_t					MethodSize = 0;
	void*					Init = NULL;
	SIMULATION_DOMAIN*		SimulationDomain = NULL;
	METHOD*					Method = NULL;
	SIMULATION*				Simulation = NULL;
	char					DomainFileName[128] = {0};
	short					NumMethods = 0;
	short					NumSims = 0;
	MSM_PARAMETERS			Parameters;
	MSM_OPTIONS				Options;

	short					pa = 4;
	short					pb = 10;
	short					plen = (pb-pa)/2 + 1;
	short					p = 0;

	short					axa = 5;
	short					axb = 20;
	short					axlen = axb-axa+1;
	short					ax = 0;

	double					a = 0.0;
	double					U = 0.0;
	double					Up = 0.0;

	//	Setup simulation domain (water sphere)
	Simulator->NumDomains = 1;
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

	//	Set up methods and simulations ("simulation" is application of method to simulation domain)
	Simulator->NumMethods = plen*axlen + 1;
	Simulator->NumSimulations = Simulator->NumDomains * Simulator->NumMethods;
	Simulator->Methods = (METHOD**) dynmem(Simulator->NumSimulations*sizeof(METHOD*));
	Simulator->Simulations = (SIMULATION**) dynmem(Simulator->NumSimulations*sizeof(SIMULATION*));
	NumMethods = 0;
	NumSims = 0;

	//	Set up NAIVE simulation
printf("\t\tSETTING UP NAIVE SIMULATION\n");
	Method = (METHOD*) dynmem(sizeof(NAIVE));
	method_initialize((void*)Method, (void*)&naive_initialize, NumMethods, NULL, NULL);
	Simulator->Methods[NumMethods] = Method;

	//	None of the previous domains have the same size as the current one
	Simulator->Methods[NumMethods]->preprocess(Simulator->Methods[NumMethods], Simulator->Domains[0]->Radius);
	NumMethods++;

	//	Create simulation (these could happen in parallel with OpenMP)
	Simulation = (SIMULATION*) dynmem(sizeof(SIMULATION));
	simulation_initialize(Simulation, Simulator->Domains[0], Simulator->Methods[NumMethods-1], NumSims, 1);	//FIXME <-- 1 is TimeSteps
	Simulator->Simulations[NumSims] = Simulation;
	NumSims++;

	//	Initialize MSM parameters
	Parameters.h = 2.5;
	Parameters.mu = 25;
	Parameters.D = 0.0;	//	Not known until preprocess/evaluate
	Parameters.L = 2;		//	# of grids

	//	Initialize MSM options
	Options.ComputeExclusions = 1;
	Options.ComputeLongRange = 1;
	Options.ComputeShortRange = 1;
	Options.IsN = 1;
	Options.IsNLogN = 0;
	Options.GridType = 0;

	//	Set up MSM simulation(s)
	for (p = pa; p <= pb; p+=2)
	{
		for (ax = axa; ax <= axb; ax++)
		{
			//	Initialize MSM parameters
			Parameters.a = ax*1.0;
			Parameters.alpha = Parameters.a / Parameters.h;
			Parameters.p = p;
			Parameters.k = p;
printf("\t\tSETTING UP MSM SIMULATION WITH p = %02hd AND a = %f\n", Parameters.p, Parameters.a);

			Method = (METHOD*) dynmem(sizeof(MSM));
			method_initialize((void*)Method, (void*)&msm_initialize, NumMethods, (void*)&Parameters, (void*)&Options);
			Simulator->Methods[NumMethods] = Method;

			//	None of the previous domains have the same size as the current one
			Simulator->Methods[NumMethods]->preprocess(Simulator->Methods[NumMethods], Simulator->Domains[0]->Radius);
			NumMethods++;

			//	Create simulation (these could happen in parallel with OpenMP)
			Simulation = (SIMULATION*) dynmem(sizeof(SIMULATION));
			simulation_initialize(Simulation, Simulator->Domains[0], Simulator->Methods[NumMethods-1], NumSims, 1);	//FIXME <-- 1 is TimeSteps
			Simulator->Simulations[NumSims] = Simulation;
			NumSims++;
		}
printf("\n");
	}

	//	Run simulations
	Simulator->NumSimulations = NumSims;
	Simulator->NumMethods = NumMethods;
	simulator_run_simulations(Simulator);

	//	Write output file
	NumSims = 0;
	U = Simulator->Simulations[NumSims]->Domain->Particles->U;
/*
	max_force = 0.0;
	for (i = 0; i < Simulator->Simulations[NumSims]->Domain->Particles->N; i++)
	{
		f = Simulator->Simulations[NumSims]->Domain->Particles->f[i];
		m = Simulator->Simulations[NumSims]->Domain->Particles->m[i];
		norm_f = sqrt(f[0]*f[0] + f[1]*f[1]+ f[2]*f[2])/m;
		if (norm_f > max_force)
			max_force = norm_f;
	}
*/
	printf("NAIVE:\t\t\t\t%f\n\n", U);
	NumSims++;

	for (p = pa; p <= pb; p+=2)
	{
		for (ax = axa; ax <= axb; ax++)
		{
			Up = Simulator->Simulations[NumSims]->Domain->Particles->U;
			printf("MSM(p=%02hd,a=%f):\t\t%f\t%e\n", p, ax*1.0, Up, fabs((U-Up)/U));
			NumSims++;
		}
		printf("\n");
	}
}

void simulator_examine_results(SIMULATOR* Simulator)
{
	double		abso = 0.0;
	double		rela = 0.0;

	long		i = 0;
	long		j = 0;
	FILE		*fp = NULL;
	char		FileName[256];

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
/*
	//	Write Energy scalar and Force vector(s) to file
	for (i = 0; i < Simulator->NumSimulations; i++)
	{
		sprintf(FileName, "%s_Sim%03ld_Results.txt", Simulator->Simulations[i]->Domain->Name, i);
		printf("Writing file: %s\n", FileName);

		//	Open file to write results
		if ((fp = fopen(FileName, "w")) != NULL)
		{
			printf("%ld\n", Simulator->Simulations[i]->Domain->Particles->N);
			fprintf(fp, "%ld\n", Simulator->Simulations[i]->Domain->Particles->N);
			printf("%.16f\n", Simulator->Simulations[i]->Domain->Particles->U);
			fprintf(fp, "%.16f\n", Simulator->Simulations[i]->Domain->Particles->U);
			for (j = 0; j < Simulator->Simulations[i]->Domain->Particles->N; j++)
			{
				printf("%.16f %.16f %.16f\n", Simulator->Simulations[i]->Domain->Particles->f[j][0],
											  Simulator->Simulations[i]->Domain->Particles->f[j][1],
											  Simulator->Simulations[i]->Domain->Particles->f[j][2]);
				fprintf(fp, "%.16f %.16f %.16f\n", Simulator->Simulations[i]->Domain->Particles->f[j][0],
												   Simulator->Simulations[i]->Domain->Particles->f[j][1],
												   Simulator->Simulations[i]->Domain->Particles->f[j][2]);
			}

			//	Close file
			fclose(fp);
		}
	}
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
	printf("ENERGY: abs O(nlogn): %e, rel O(nlogn): %e\n", abso, rela);
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
	printf("FORCE: MAX abs O(nlong): %e, rel O(nlogn): %e\n", abso, -1.0);
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