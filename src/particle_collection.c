//-------|---------|---------|---------|---------|---------|---------|---------|
/*
particle_collection.c - 
*/

#include "particle_collection.h"

//	EXTERNAL Methods
void particle_collection_initialize(PARTICLE_COLLECTION* Pc, long N, double UnitConverter)
{
	//	A particle collection is a list of particles and any relevant properties (mass, velocity, etc)
	assert(Pc != NULL);

	Pc->N = N;
	Pc->U = 0.0;
	Pc->r = (PARTICLE*) dynmem(Pc->N*sizeof(PARTICLE));
	Pc->f = (double**) dynarr_d(Pc->N, 3);
	Pc->m = (double*) dynvec(Pc->N, sizeof(double));
	Pc->v = (double**) dynarr_d(Pc->N, 3);
	Pc->UnitConverter = UnitConverter;
}

void particle_collection_update(PARTICLE_COLLECTION* Pc, PARTICLE* MinimumCoordinate, PARTICLE* MaximumCoordinate)
{
	long		i = 0;

	//	This method will update the particle positions at each time step
	//		->	Leap-frog integration

	MinimumCoordinate->x = DBL_MAX;
	MinimumCoordinate->y = DBL_MAX;
	MinimumCoordinate->z = DBL_MAX;
	MinimumCoordinate->q = DBL_MAX;
	MaximumCoordinate->x = -DBL_MAX;
	MaximumCoordinate->y = -DBL_MAX;
	MaximumCoordinate->z = -DBL_MAX;
	MaximumCoordinate->q = -DBL_MAX;

	//	Loop through all of the particles and update their positions
	for (i = 0; i < Pc->N; i++)
	{
		//	New position = old position + momemtum + forces + ...
		//		... FIXME

		//	Check particle i to see if it has a minimum coordinate
		MinimumCoordinate->x = MIN(MinimumCoordinate->x, Pc->r[i].x);
		MinimumCoordinate->y = MIN(MinimumCoordinate->y, Pc->r[i].y);
		MinimumCoordinate->z = MIN(MinimumCoordinate->z, Pc->r[i].z);
		MinimumCoordinate->q = MIN(MinimumCoordinate->q, Pc->r[i].q);
		//	Check particle i to see if it has a maximum coordinate
		MaximumCoordinate->x = MAX(MaximumCoordinate->x, Pc->r[i].x);
		MaximumCoordinate->y = MAX(MaximumCoordinate->y, Pc->r[i].y);
		MaximumCoordinate->z = MAX(MaximumCoordinate->z, Pc->r[i].z);
		MaximumCoordinate->q = MAX(MaximumCoordinate->q, Pc->r[i].q);
	}
}

void particle_collection_copy(PARTICLE_COLLECTION* DstParticles, PARTICLE_COLLECTION* SrcParticles)
{
	long i = 0;

	assert(SrcParticles != NULL);
	assert(DstParticles != NULL);

	//	Create template for particle_collection
	particle_collection_initialize(DstParticles, SrcParticles->N, SrcParticles->UnitConverter);

	//	Energy, U, is initialized to 0.0. Copy over the correct value from SrcParticles
	DstParticles->U = SrcParticles->U;

	//	Fill in all of the particle info from source to destination
	for (i = 0; i < DstParticles->N; i++)
	{
		DstParticles->r[i].x = SrcParticles->r[i].x;
		DstParticles->r[i].y = SrcParticles->r[i].y;
		DstParticles->r[i].z = SrcParticles->r[i].z;
		DstParticles->r[i].q = SrcParticles->r[i].q;

		DstParticles->f[i][0] = SrcParticles->f[i][0];
		DstParticles->f[i][1] = SrcParticles->f[i][1];
		DstParticles->f[i][2] = SrcParticles->f[i][2];

		DstParticles->m[i] = SrcParticles->m[i];

		DstParticles->v[i][0] = SrcParticles->v[i][0];
		DstParticles->v[i][1] = SrcParticles->v[i][1];
		DstParticles->v[i][2] = SrcParticles->v[i][2];
/*
		printf("Particle %+02ld: (%+f,%+f,%+f), q=%+f,  m=%+f, v=(%+f,%+f,%+f)\n",
			i, DstParticles->r[i].x, DstParticles->r[i].y, DstParticles->r[i].z, DstParticles->r[i].q,
			DstParticles->m[i],
			DstParticles->v[i][0], DstParticles->v[i][1], DstParticles->v[i][2]);
*/
	}
}

void particle_collection_uninitialize(PARTICLE_COLLECTION* Pc)
{
	//	Free dynamically allocated memory
printf("before dynfree(Pc->v[0])\n");
	dynfree(Pc->v[0]);
printf("before dynfree(Pc->v)\n");
	dynfree(Pc->v);
printf("before dynfree(Pc->m)\n");
	dynfree(Pc->m);
printf("before dynfree(Pc->f[0])\n");
	dynfree(Pc->f[0]);
printf("before dynfree(Pc->f)\n");
	dynfree(Pc->f);
printf("before dynfree(Pc->r)\n");
	dynfree(Pc->r);
}

//	End of file
