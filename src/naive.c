//-------|---------|---------|---------|---------|---------|---------|---------|
/*
naive.c -
*/

#include "naive.h"

//	EXTERNAL Methods
void naive_initialize(void* Method)
{
	NAIVE*		Naive = (NAIVE*) Method;

	assert(Naive != NULL);
	printf("Initializing NAIVE!\n");

	//	Initialize COMMON members
	Naive->cmn.Size = sizeof(NAIVE);

	//	Initialize COMMON function pointers
	Naive->cmn.copy = &naive_copy;
	Naive->cmn.preprocess = &naive_preprocess;
	Naive->cmn.evaluate = &naive_evaluate;
	Naive->cmn.uninitialize = &naive_uninitialize;
}

void naive_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	--> METHOD is copied in method_copy()
}

void naive_preprocess(void* Method, double DomainRadius)
{
	NAIVE*		Naive = (NAIVE*) Method;
	assert(Naive != NULL);
	printf("NAIVE Preprocessing!\n");
}

void naive_evaluate(void* Method, long N, PARTICLE* r, double* U, double** f)
{
	NAIVE*		Naive = (NAIVE*) Method;
	long		i = 0;
	long		j = 0;
	double		d = 0.0;
	double		dx = 0.0;
	double		dy = 0.0;
	double		dz = 0.0;
	double		dfx = 0.0;
	double		dfy = 0.0;
	double		dfz = 0.0;

	assert(Naive != NULL);
	printf("NAIVE Evaluation! %lu particles\n", N);

	//	Perform the naive O(N^2) calculation
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			//	Compute Euclidean distance between two particles
			dx = r[i].x-r[j].x;
			dy = r[i].y-r[j].y;
			dz = r[i].z-r[j].z;
			d = sqrt(dx*dx + dy*dy + dz*dz);

			//	Compute contribution to the energy
			*U += (r[i].q*r[j].q/d);

			//	Compute contribution to the forces
			dfx = (-dx/d)*r[i].q*r[j].q/(d*d);
			dfy = (-dy/d)*r[i].q*r[j].q/(d*d);
			dfz = (-dz/d)*r[i].q*r[j].q/(d*d);

			//	Apply force to particle i
			f[i][0] -= dfx;
			f[i][1] -= dfy;
			f[i][2] -= dfz;

			//	Apply force to particle j
			f[j][0] += dfx;
			f[j][1] += dfy;
			f[j][2] += dfz;
		}
	}
}

void naive_uninitialize(void* Method)
{
	NAIVE*		Naive = (NAIVE*) Method;
	assert(Naive != NULL);
	printf("Un-initializing NAIVE!\n");
}

//	INTERNAL Methods

//	End of file
