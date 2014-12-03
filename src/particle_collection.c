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
	Pc->r = (PARTICLE*) dynmem(Pc->N*sizeof(PARTICLE));
	Pc->m = (double*) dynvec(Pc->N, sizeof(double));
	Pc->v = (double**) dynarr_d(Pc->N, 3);
	Pc->UnitConverter = UnitConverter;
}

void particle_collection_move(PARTICLE_COLLECTION* Pc)
{
	//	This method will update the particle positions at each time step
}

void particle_collection_copy(PARTICLE_COLLECTION* SrcParticles, PARTICLE_COLLECTION* DstParticles)
{
	long i = 0;

	assert(SrcParticles != NULL);
	assert(DstParticles != NULL);

	particle_collection_initialize(DstParticles, SrcParticles->N, SrcParticles->UnitConverter);

	//DstParticles->N = SrcParticles->N;
	//DstParticles->UnitConverter = SrcParticles->UnitConverter;
	//memcpy(&DstParticles->r, &SrcParticles->r, DstParticles->N*sizeof(PARTICLE));
	//memcpy(&DstParticles->m, &SrcParticles->m, DstParticles->N*sizeof(double));
	//memcpy(&DstParticles->v, &SrcParticles->v, DstParticles->N*3*sizeof(double)); // does this make sense?
	for (i = 0; i < DstParticles->N; i++)
	{
		DstParticles->r[i].x = SrcParticles->r[i].x;
		DstParticles->r[i].y = SrcParticles->r[i].y;
		DstParticles->r[i].z = SrcParticles->r[i].z;
		DstParticles->r[i].q = SrcParticles->r[i].q;

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
printf("before dynfree(Pc->r)\n");
	dynfree(Pc->r);
}

//	End of file
