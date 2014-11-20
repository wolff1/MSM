//-------|---------|---------|---------|---------|---------|---------|---------|
/*
domain.c - 
*/

#include "simulation_domain.h"

//	EXTERNAL Methods
void simulation_domain_initialize(SIMULATION_DOMAIN** Domain, short Id, char* FileName)
{
	//	A domain is a system of particles, and any other relevant information
	assert(*Domain == NULL);

	(*Domain) = (SIMULATION_DOMAIN*) dynmem(sizeof(SIMULATION_DOMAIN));

	(*Domain)->Id = Id;
	strncpy((*Domain)->Name, FileName, MAXLEN_DOMAIN_NAME_STR);
	(*Domain)->MinimumCoordinates.x = DBL_MAX;
	(*Domain)->MinimumCoordinates.y = DBL_MAX;
	(*Domain)->MinimumCoordinates.z = DBL_MAX;
	(*Domain)->MaximumCoordinates.x = -DBL_MAX;
	(*Domain)->MaximumCoordinates.y = -DBL_MAX;
	(*Domain)->MaximumCoordinates.z = -DBL_MAX;

	simulation_domain_input_particles(*Domain);

//	printf("MIN = (%f,%f,%f)\n", (*Domain)->MinimumCoordinates.x,(*Domain)->MinimumCoordinates.y,(*Domain)->MinimumCoordinates.z);
//	printf("CTR = (%f,%f,%f)\n", (*Domain)->CenterCoordinates.x,(*Domain)->CenterCoordinates.y,(*Domain)->CenterCoordinates.z);
//	printf("MAX = (%f,%f,%f)\n", (*Domain)->MaximumCoordinates.x,(*Domain)->MaximumCoordinates.y,(*Domain)->MaximumCoordinates.z);
//	printf("Radius = %f, N=%ld\n", (*Domain)->Radius, (*Domain)->Particles->N);
}

void simulation_domain_compute_forces(SIMULATION_DOMAIN* Domain)
{
	//	Pass domain info to method and let it calculate the forces and electrostatic energy
}

void simulation_domain_uninitialize(SIMULATION_DOMAIN* Domain)
{
	//	Free dynamically allocated memory
	particle_collection_uninitialize(Domain->Particles);
	dynfree(Domain);
}

//	INTERNAL Methods
void simulation_domain_input_particles(SIMULATION_DOMAIN* Domain)
{
	long						i = 0;
	long						Idx = 0;
	FILE*						Fp = NULL;
	PARTICLE_COLLECTION*		Pc = NULL;
	double						UnitConverter = 1.0;
	char						FileName[132] = {0};
	char						TempBuffer[128] = {0};
	FILE_FORMAT_PSF				LinePsf;
	FILE_FORMAT_PDB				LinePdb;
	long						N = 0;

	if (1)
	{	//	kcal/mol
		UnitConverter = sqrt(UNITS_COULOMB);
	}

	if (Domain->Name != NULL)
	{
		//	Read PSF file, set up Particles
		sprintf(FileName, "%s.psf", Domain->Name);
		Fp = fopen(FileName, "r");

		while (1)
		{
			memset(&LinePsf, 0, sizeof(FILE_FORMAT_PSF));

			//	Read lines until we reach !NATOM or end of file
			if (fgets((char*)&LinePsf, sizeof(FILE_FORMAT_PSF), Fp) == NULL)
			{
				break;
			}

			if (!strncmp(LinePsf.Header.Type, "!NATOM", 6))
			{
				strncpy(TempBuffer, LinePsf.Header.NumRecs, 8);
				N = strtol(TempBuffer, NULL, 10);
				break;
			}
		}

		//	Now that we know how many particles, set up particle_collection
		particle_collection_initialize(&Pc, N, UnitConverter);

		//	Read rest of file(s)
		for (i = 0; i < Pc->N; i++)
		{
			memset(&LinePsf, 0, sizeof(FILE_FORMAT_PSF));
			memset(&TempBuffer, 0, 128);

			if (fgets((char*)&LinePsf, sizeof(FILE_FORMAT_PSF), Fp) == NULL)
			{
				break;
			}

			strncpy(TempBuffer, LinePsf.Atom.UID, 8);
			Idx = strtol(TempBuffer, NULL, 10) - 1;

			if (Idx < Pc->N)
			{
				strncpy(TempBuffer, LinePsf.Atom.Charge, 16);
				Pc->r[Idx].q = strtod(TempBuffer, NULL);

				strncpy(TempBuffer, LinePsf.Atom.Mass, 7);
				Pc->m[Idx] = strtod(TempBuffer, NULL);
			}
		}

		fclose(Fp);

		//	Read PDB file, finish setting up Particles
		sprintf(FileName, "%s.pdb", Domain->Name);
		Fp = fopen(FileName, "r");

		while (1)
		{
			memset(&LinePdb, 0, sizeof(FILE_FORMAT_PDB));
			memset(&TempBuffer, 0, 128);

			if(fgets((char*)&LinePdb, sizeof(FILE_FORMAT_PDB), Fp) == NULL)
			{
				break;
			}

			//	Ignore line if it is NOT an ATOM record.
			if (strncmp(LinePdb.Atom.RecName, "ATOM", 4))
			{
				continue;
			}

			strncpy(TempBuffer, LinePdb.Atom.UID, 5);
			Idx = strtol(TempBuffer, NULL, 10) - 1;

			if (Idx < Pc->N)
			{
				strncpy(TempBuffer, LinePdb.Atom.X, 8);
				Pc->r[Idx].x = strtod(TempBuffer, NULL);

				strncpy(TempBuffer, LinePdb.Atom.Y, 8);
				Pc->r[Idx].y = strtod(TempBuffer, NULL);

				strncpy(TempBuffer, LinePdb.Atom.Z, 8);
				Pc->r[Idx].z = strtod(TempBuffer, NULL);

				Domain->MinimumCoordinates.x = MIN(Domain->MinimumCoordinates.x, Pc->r[Idx].x);
				Domain->MinimumCoordinates.y = MIN(Domain->MinimumCoordinates.y, Pc->r[Idx].y);
				Domain->MinimumCoordinates.z = MIN(Domain->MinimumCoordinates.z, Pc->r[Idx].z);

				Domain->MaximumCoordinates.x = MAX(Domain->MaximumCoordinates.x, Pc->r[Idx].x);
				Domain->MaximumCoordinates.y = MAX(Domain->MaximumCoordinates.y, Pc->r[Idx].y);
				Domain->MaximumCoordinates.z = MAX(Domain->MaximumCoordinates.z, Pc->r[Idx].z);
			}
		}

		fclose(Fp);
	}
	else
	{
		//	Generate data, set up Particles
	}

	//	Find lengths of sides of simulation domain
	Domain->CenterCoordinates.x = Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x;
	Domain->CenterCoordinates.y = Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y;
	Domain->CenterCoordinates.z = Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z;

	//	Define Radius to be 2-norm of vector between the min and max corners
	Domain->Radius =  Domain->CenterCoordinates.x * Domain->CenterCoordinates.x;
	Domain->Radius += Domain->CenterCoordinates.y * Domain->CenterCoordinates.y;
	Domain->Radius += Domain->CenterCoordinates.z * Domain->CenterCoordinates.z;
	Domain->Radius = 0.5 * sqrt(Domain->Radius);

	//	Lastly, set the center to be the minimum coordinate plus half of the side length
	Domain->CenterCoordinates.x = Domain->MinimumCoordinates.x + (0.5 * Domain->CenterCoordinates.x);
	Domain->CenterCoordinates.y = Domain->MinimumCoordinates.y + (0.5 * Domain->CenterCoordinates.y);
	Domain->CenterCoordinates.z = Domain->MinimumCoordinates.z + (0.5 * Domain->CenterCoordinates.z);

	Domain->Particles = Pc;
}

//	End of file