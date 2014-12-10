//-------|---------|---------|---------|---------|---------|---------|---------|
/*
msm.c - 
*/

#include "msm.h"

//	EXTERNAL Methods
void msm_initialize(void* Method)
{
	size_t		Size = 0;
	void*		Init = NULL;
	void*		Ptr = NULL;
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	printf("Initializing MSM!\n");

	//	Initialize COMMON members
	Msm->cmn.Size = sizeof(MSM);

	//	Initialize COMMON function pointers
	Msm->cmn.copy = &msm_copy;
	Msm->cmn.preprocess = &msm_preprocess;
	Msm->cmn.evaluate = &msm_evaluate;
	Msm->cmn.uninitialize = &msm_uninitialize;

	//	Initialize MSM parameters
	Msm->prm.a = 12.5;
	Msm->prm.h = 2.5;
	Msm->prm.alpha = Msm->prm.a / Msm->prm.h;
	Msm->prm.p = 4;
	Msm->prm.k = 4;
	Msm->prm.mu = 10;
	Msm->prm.D = 0.0;	//	Not known until preprocess/evaluate

	//	Initialize MSM options
	Msm->opt.ComputeExclusions = 1;
	Msm->opt.ComputeLongRange = 1;
	Msm->opt.ComputeShortRange = 1;
	Msm->opt.IsN = 1;
	Msm->opt.IsNLogN = 1;

	//	Initialize INTERPOLANT
	Ptr = NULL;
	Msm->itp = NULL;
	if (1)
	{	//	B_SPLINE interpolant
		Size = sizeof(B_SPLINE);
		Init = &b_spline_initialize;
	}
	else
	{	//	C1_SPLINE interpolant
		Size = sizeof(C1_SPLINE);
		Init = &c1_spline_initialize;
	}
	Ptr = (INTERPOLANT*) dynmem(Size);
	interpolant_initialize(Ptr, Init, &Msm->prm);
	Msm->itp = (INTERPOLANT*) Ptr;

	//	Initialize SOFTENER
	Ptr = NULL;
	Msm->sft = NULL;
	if (1)
	{	//	EVEN_POWERS softening (aka "Taylor")
		Size = sizeof(EVEN_POWERS);
		Init = &even_powers_initialize;
	}
	Ptr = (SOFTENER*) dynmem(Size);
	softener_initialize(Ptr, Init, Msm->prm.k);
	Msm->sft = (SOFTENER*) Ptr;
}

void msm_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	--> METHOD is copied in method_copy()

	//	Copy MSM_PARAMETERS
	memcpy(&((MSM*)Dst)->prm, &((MSM*)Src)->prm, sizeof(MSM_PARAMETERS));

	//	Copy MSM_OPTIONS
	memcpy(&((MSM*)Dst)->opt, &((MSM*)Src)->opt, sizeof(MSM_OPTIONS));

	//	Copy INTERPOLANT
	((MSM*)Dst)->itp = (INTERPOLANT*) dynmem(((MSM*)Src)->itp->Size);

	interpolant_copy(((MSM*)Dst)->itp, ((MSM*)Src)->itp);

	//	Copy SOFTENER
	((MSM*)Dst)->sft = (SOFTENER*) dynmem(((MSM*)Src)->sft->Size);
	softener_copy(((MSM*)Dst)->sft, ((MSM*)Src)->sft);
}

void msm_preprocess(void* Method, double DomainRadius)
{
	MSM*		Msm = (MSM*) Method;
	assert(Msm != NULL);
	printf("MSM Preprocessing for size <%f>!\n", DomainRadius);

	//	Compute Interpolant coefficients
	Msm->prm.D = DomainRadius;
	(*Msm->itp->compute_g2g)(Msm->itp, Msm->sft, &Msm->prm);
	(*Msm->itp->compute_tg2g)(Msm->itp, Msm->sft, &Msm->prm);
}

void msm_evaluate(void* Method, SIMULATION_DOMAIN* Domain)
{
	long		N = 0;
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	assert(Domain != NULL);

	N = Domain->Particles->N;
	printf("MSM Evaluation! %lu particles\n", N);

	//	Initialize output variables U and f
	Domain->Particles->U = 0.0;
	memset(Domain->Particles->f[0], 0, sizeof(double)*N*3);	//	FIXME <-- double check this

	if (Msm->opt.ComputeShortRange)
	{
		msm_short_range(Msm, Domain);
	}

	if (Msm->opt.ComputeLongRange)
	{
		if (Msm->opt.IsN)
		{
			msm_anterpolate(Msm);
			msm_restrict(Msm);
			msm_direct(Msm);
			msm_direct_top(Msm);
			msm_prolongate(Msm);
			msm_interpolate(Msm);
		}

		if (Msm->opt.IsNLogN)
		{
			msm_anterpolate(Msm);
			msm_direct_top(Msm);
			msm_interpolate(Msm);
		}

		if (Msm->opt.ComputeExclusions)
		{
			msm_exclude(Msm);
		}
	}
}

void msm_uninitialize(void* Method)
{
	MSM*		Msm = (MSM*) Method;
	assert(Msm != NULL);
	printf("Un-initializing MSM!\n");

	//	Uninitialize SOFTENER
	(*Msm->sft->uninitialize)(Msm->sft);
	dynfree(Msm->sft);

	//	Uninitialize INTERPOLANT
	(*Msm->itp->uninitialize)(Msm->itp);
	dynfree(Msm->itp);
}

//	INTERNAL Methods
long IDX(long x, long y, long z, long YLength, long ZLength)
{
	return x*YLength*ZLength + y*ZLength + z;
}

void msm_short_range(MSM* Msm, SIMULATION_DOMAIN* Domain)
{
	long		i = 0;
	long		j = 0;
	long		k = 0;
	long		n = 0;
	long		N = Domain->Particles->N;
	PARTICLE*	r = Domain->Particles->r;
	double*		U = &Domain->Particles->U;
	double**	f = Domain->Particles->f;
	double		a = Msm->prm.a;
	long		XBinCount = 0;
	long		YBinCount = 0;
	long		ZBinCount = 0;
	long*		First = NULL;
	long*		Next = NULL;

	printf("\tMSM short range computation!\n");

	//	Set up number of bins in each dimension (NOTE: MATLAB code has floor() + 1)
	XBinCount = (long) ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / a);
	YBinCount = (long) ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / a);
	ZBinCount = (long) ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / a);

	//	Dynamically allocate memory for lists of bins
	First = (long*) dynvec(XBinCount*YBinCount*ZBinCount, sizeof(long));
	Next = (long*) dynvec(N, sizeof(long));

	//	Set-up lists of bins
	for (n = 0; n < XBinCount*YBinCount*ZBinCount; n++)
	{
		First[n] = -1;
	}

	//	Assign particles to bins
	for (n = 0; n < N; n++)
	{
		i = (long) floor((r[n].x-Domain->MinimumCoordinates.x) / a);
		j = (long) floor((r[n].y-Domain->MinimumCoordinates.y) / a);
		k = (long) floor((r[n].z-Domain->MinimumCoordinates.z) / a);
		Next[n] = First[IDX(i,j,k,YBinCount,ZBinCount)];
		First[IDX(i,j,k,YBinCount,ZBinCount)] = n;
	}

	//	Loop over bins
	for (i = 0; i < XBinCount; i++)
	{
		for (j = 0; j < YBinCount; j++)
		{
			for (k = 0; k < ZBinCount; k++)
			{
//				printf("%03ld <-- (%02ld,%02ld,%02ld)\n", IDX(i,j,k,YBinCount,ZBinCount), i,j,k);
				if (First[IDX(i,j,k,YBinCount,ZBinCount)] != -1)
				{
					printf("Bin: (%02ld,%02ld,%02ld):\t", i,j,k);
					n = First[IDX(i,j,k,YBinCount,ZBinCount)];
					do
					{
						printf("%02ld, ", n);
						n = Next[n];
					} while (n != -1);
					printf("\n");
				}
			}
		}
	}

	//	Free dynamically allocated memory
	dynfree(First);
	dynfree(Next);
}

void msm_anterpolate(MSM* Msm)
{
}

void msm_restrict(MSM* Msm)
{
}

void msm_direct(MSM* Msm)
{
}

void msm_direct_top(MSM* Msm)
{
}

void msm_prolongate(MSM* Msm)
{
}

void msm_interpolate(MSM* Msm)
{
}

void msm_exclude(MSM* Msm)
{
}

//	End of file
