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
	Msm->prm.a = 1.5;
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
	Msm->opt.IsN = 0;
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
	memset(Domain->Particles->f[0], 0, sizeof(double)*N*3);

	if (Msm->opt.ComputeShortRange)
	{
		msm_short_range(Msm, Domain);
	}

	if (Msm->opt.ComputeLongRange)
	{
		if (Msm->opt.IsN)
		{
			msm_anterpolate(Msm, Domain, 0);
			msm_restrict(Msm);
			msm_direct(Msm);
			msm_direct_top(Msm);
			msm_prolongate(Msm);
			msm_interpolate(Msm);
		}

		if (Msm->opt.IsNLogN)
		{
			for (N = 0; N < 3; N++)
			{	//	FIXME
				msm_anterpolate(Msm, Domain, (short)N);
				msm_direct_top(Msm);
				msm_interpolate(Msm);
			}
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
	long		Particle1 = 0;
	long		Particle2 = 0;
	long		MaxInteractionCount = 0;
	long*		First = NULL;
	long*		Next = NULL;
	long*		ParticlesPerBin = NULL;
	long		ParticlesInBin = 0;
	short		EdgeCaseCount = 0;

	printf("\tMSM short range computation!\n");

	//	Set up number of bins in each dimension (NOTE: MATLAB code has floor() + 1)
	XBinCount = (long) ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / a);
	YBinCount = (long) ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / a);
	ZBinCount = (long) ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / a);
printf("\t\t%ldx%ldx%ld bins\n", XBinCount, YBinCount, ZBinCount);
	//	Dynamically allocate memory for lists of bins
	First = (long*) dynvec(XBinCount*YBinCount*ZBinCount, sizeof(long));
	Next = (long*) dynvec(N, sizeof(long));
	ParticlesPerBin = (long*) dynvec(XBinCount*YBinCount*ZBinCount, sizeof(long));

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
		Next[n] = First[IDX(i,j,k,XBinCount,YBinCount)];
		First[IDX(i,j,k,XBinCount,YBinCount)] = n;
		ParticlesPerBin[IDX(i,j,k,XBinCount,YBinCount)]++;
	}
/*
	//	Loop over bins
	for (k = 0; k < ZBinCount; k++)
	{
		for (j = 0; j < YBinCount; j++)
		{
			for (i = 0; i < XBinCount; i++)
			{
				printf("\t\t%03ld <-- (%02ld,%02ld,%02ld), Cnt=%ld\n", IDX(i,j,k,XBinCount,YBinCount), i,j,k, ParticlesPerBin[IDX(i,j,k,XBinCount,YBinCount)]);
				if (First[IDX(i,j,k,XBinCount,YBinCount)] != -1)
				{
//					printf("Bin: (%02ld,%02ld,%02ld):\t", i,j,k);
					n = First[IDX(i,j,k,XBinCount,YBinCount)];
					do
					{
						printf("%02ld, ", n);
						n = Next[n];
					} while (n != -1);
					printf("\n");
				}
				printf("\n");
			}
		}
	}
*/

	for (k = 1; k < ZBinCount; k++)
	{
		for (j = 1; j < YBinCount; j++)
		{
			for (i = 1; i < XBinCount; i++)
			{
				//	"next" bin is (i-1, j-1, k-1) -> i > 0, j > 0, k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i-1,j-1,k-1);
			}

			for (i = 0; i < XBinCount; i++)
			{
				//	"next" bin is (i, j-1, k-1) ->        j > 0, k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i,j-1,k-1);
			}

			for (i = 0; i < XBinCount-1; i++)
			{
				//	"next" bin is (i+1, j-1, k-1) ->	i+1<X, j > 0, k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i+1,j-1,k-1);
			}
		}

		for (j = 0; j < YBinCount; j++)
		{
			for (i = 1; i < XBinCount; i++)
			{
				//	"next" bin is (i-1, j, k-1) ->	i > 0,        k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i-1,j,k-1);
			}

			for (i = 0; i < XBinCount; i++)
			{
				//	"next" bin is (i, j, k-1) ->	              k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i,j,k-1);
			}

			for (i = 0; i < XBinCount-1; i++)
			{
				//	"next" bin is (i+1, j, k-1) ->	i+1<X,        k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i+1,j,k-1);
			}
		}

		for (j = 0; j < YBinCount-1; j++)
		{
			for (i = 1; i < XBinCount; i++)
			{
				//	"next" bin is (i-1, j+1, k-1) ->	i > 0, j+1<Y, k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i-1,j+1,k-1);
			}

			for (i = 0; i < XBinCount; i++)
			{
				//	"next" bin is (i, j+1, k-1) ->	       j+1<Y, k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i,j+1,k-1);
			}

			for (i = 0; i < XBinCount-1; i++)
			{
				//	"next" bin is (i+1, j+1, k-1) ->	i+1<X, j+1<Y, k > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i+1,j+1,k-1);
			}
		}
	}

	for (k = 0; k < ZBinCount; k++)
	{
		for (j = 1; j < YBinCount; j++)
		{
			for (i = 1; i < XBinCount; i++)
			{
				//	"next" bin is (i-1, j-1, k) ->	i > 0, j > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i-1,j-1,k);
			}

			for (i = 0; i < XBinCount; i++)
			{
				//	"next" bin is (i, j-1, k) ->	       j > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i,j-1,k);
			}

			for (i = 0; i < XBinCount-1; i++)
			{
				//	"next" bin is (i+1, j-1, k) ->	i+1<X, j > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i+1,j-1,k);
			}
		}

		for (j = 0; j < YBinCount; j++)
		{
			for (i = 1; i < XBinCount; i++)
			{
				//	"next" bin is (i-1, j, k) ->	i > 0
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i,j,k, i-1,j,k);
			}

			for (i = 0; i < XBinCount; i++)
			{
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, XBinCount, YBinCount, i, j, k);
			}
		}
	}

	//	Free dynamically allocated memory
	dynfree(First);
	dynfree(Next);
	dynfree(ParticlesPerBin);
}

void msm_anterpolate(MSM* Msm, SIMULATION_DOMAIN* Domain, short Level)
{
	long		n = 0;
	long		Nx = 0;
	long		Ny = 0;
	long		Nz = 0;
	double		h = Msm->prm.h;
	short		p = Msm->prm.p;

	//	Create grid corresponding to Level (0 is finest grid)
	Nx = (long) ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / (h*pow(2.0, Level))) + p + 1;
	Ny = (long) ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / (h*pow(2.0, Level))) + p + 1;
	Nz = (long) ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / (h*pow(2.0, Level))) + p + 1;

	printf("Min(%f,%f,%f) Max(%f,%f,%f) -> %ldx%ldx%ld (Level=%hd, h=%2.2f)\n",
		Domain->MinimumCoordinates.x,Domain->MinimumCoordinates.y,Domain->MinimumCoordinates.z,
		Domain->MaximumCoordinates.x,Domain->MaximumCoordinates.y,Domain->MaximumCoordinates.z,
		Nx,Ny,Nz,Level,h*pow(2.0, Level));

	//	FIXME
	//	Think about local vs global indexing of grid
	//	Think about how consecutive grids will be aligned
	//	Think about how p/2 extra points on ends affect indexing
	//	Where to store grid(s)?
	//	What type of meta-data to store with grid(s)?

	//	Loop through all particles
	//		-> Spread particle charge onto grid in each dimension
	for (n = 0; n < Domain->Particles->N; n++)
	{
	}
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

//	INTERNAL HELPER Methods
void msm_short_range_bin_to_bin(MSM* Msm, SIMULATION_DOMAIN* Domain, long* Next, long Particle1, long Particle2, long MaxIterationCount)
{
	long		i = Particle1;
	long		j = Particle2;
	long		k = 0;
	double		a = Msm->prm.a;
	double		a2 = a*a;
	PARTICLE*	r = Domain->Particles->r;
	double*		U = &Domain->Particles->U;
	double**	f = Domain->Particles->f;
	double		dx = 0.0;
	double		dy = 0.0;
	double		dz = 0.0;
	double		d = 0.0;
	long		Idx = 0;
	double*		D = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	PARTICLE*	R = NULL;
	long*		I = NULL;
	long*		J = NULL;
	double		Magnitude = 0.0;
	double		Direction[3];

	//	Dynamically allocate memory
	D = (double*) dynvec(MaxIterationCount, sizeof(double));
	F = (double*) dynvec(MaxIterationCount, sizeof(double));
	DF = (double*) dynvec(MaxIterationCount, sizeof(double));
	R = (PARTICLE*) dynvec(MaxIterationCount, sizeof(PARTICLE));
	I = (long*) dynvec(MaxIterationCount, sizeof(long));
	J = (long*) dynvec(MaxIterationCount, sizeof(long));

	//	Particle i in bin 1
	while (i != -1)
	{
		if (Particle1 != Particle2)
		{
			j = Particle2;	//	Resets particle j to first in the list for new particle i (want this if the bins are different)
		}
		else
		{
			j = Next[i];	//	Resets particle j to first in the list for new particle i (want this if the bins are the same)
		}

		//	Particle j in bin 2
		while (j != -1)
		{
			dx = r[j].x - r[i].x;
			dy = r[j].y - r[i].y;
			dz = r[j].z - r[i].z;
			d = sqrt(dx*dx + dy*dy + dz*dz);

			if (d < a)
			{
				//	Add d/a to vector D, NOTE: need a*D[x] to get d below
				D[Idx] = d/a;

				//	Add r[j]-r[i] direction vector and q[j]*q[i] to array R
				R[Idx].x = dx;
				R[Idx].y = dy;
				R[Idx].z = dz;
				R[Idx].q = r[j].q * r[i].q;

				//	Add i to vector I and add j to vector J
				I[Idx] = i;
				J[Idx] = j;

				Idx++;
			}
			j = Next[j];
		}
		i = Next[i];
	}

	//	Compute gamma(D) and gamma'(D) in bulk
	(*Msm->sft->soften)(Msm->sft, Idx, D, F, DF);

	//	Compute energy and forces
	for (k = 0; k < Idx; k++)
	{
		//	U += q(i)*q(j)*(1/d - gamma(d/a)/a);
		(*U) += R[k].q*(1.0/(a*D[k]) - F[k]/a);

		//	f += q(i)*q(j)*(1/d^2 + gamma'(d/a)/a^2) * (r[j]-r[i])/d;
		Magnitude = R[k].q*(1.0/(a2*D[k]*D[k]) + DF[k]/a2);
		Direction[0] = R[k].x / (a*D[k]);
		Direction[1] = R[k].y / (a*D[k]);
		Direction[2] = R[k].z / (a*D[k]);

		f[I[k]][0] -= Magnitude*Direction[0];
		f[I[k]][1] -= Magnitude*Direction[1];
		f[I[k]][2] -= Magnitude*Direction[2];

		f[J[k]][0] += Magnitude*Direction[0];
		f[J[k]][1] += Magnitude*Direction[1];
		f[J[k]][2] += Magnitude*Direction[2];
	}

	//	Free dynamically allocated memory
	dynfree(D);
	dynfree(F);
	dynfree(DF);
	dynfree(R);
	dynfree(I);
	dynfree(J);
}

void msm_short_range_compute_self(MSM* Msm, SIMULATION_DOMAIN* Domain, long* First, long* Next, long* ParticlesPerBin, long XBinCount, long YBinCount, long x, long y, long z)
{
	long		i = -1;
	long		j = -1;
	long		k = 0;
	double		a = Msm->prm.a;
	double		a2 = a*a;
	PARTICLE*	r = Domain->Particles->r;
	double*		U = &Domain->Particles->U;
	double**	f = Domain->Particles->f;
	double		dx = 0.0;
	double		dy = 0.0;
	double		dz = 0.0;
	double		d = 0.0;
	long		Idx = 0;
	double*		D = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	PARTICLE*	R = NULL;
	long*		I = NULL;
	long*		J = NULL;
	double		Magnitude = 0.0;
	double		Direction[3];
	long		ParticlesInBin = 0;
	long		MaxInteractionCount = 0;

	//	What is the maximum number of interactions which we could compute?
	ParticlesInBin = ParticlesPerBin[IDX(x,y,z,XBinCount,YBinCount)];
	MaxInteractionCount = (ParticlesInBin*(ParticlesInBin+1))/2;

	//	Dynamically allocate memory
	D = (double*) dynvec(MaxInteractionCount, sizeof(double));
	F = (double*) dynvec(MaxInteractionCount, sizeof(double));
	DF = (double*) dynvec(MaxInteractionCount, sizeof(double));
	R = (PARTICLE*) dynvec(MaxInteractionCount, sizeof(PARTICLE));
	I = (long*) dynvec(MaxInteractionCount, sizeof(long));
	J = (long*) dynvec(MaxInteractionCount, sizeof(long));

	//	Particle i in bin 1
	i = First[IDX(x,y,z,XBinCount,YBinCount)];
	while (i != -1)
	{
		j = Next[i];	//	Resets particle j to first in the list for new particle i

		//	Particle j in bin 2
		while (j != -1)
		{
			dx = r[j].x - r[i].x;
			dy = r[j].y - r[i].y;
			dz = r[j].z - r[i].z;
			d = sqrt(dx*dx + dy*dy + dz*dz);

			if (d < a)
			{
				//	Add d/a to vector D, NOTE: need a*D[x] to get d below
				D[Idx] = d/a;

				//	Add r[j]-r[i] direction vector and q[j]*q[i] to array R
				R[Idx].x = dx;
				R[Idx].y = dy;
				R[Idx].z = dz;
				R[Idx].q = r[j].q * r[i].q;

				//	Add i to vector I and add j to vector J
				I[Idx] = i;
				J[Idx] = j;

				Idx++;
			}
			j = Next[j];
		}
		i = Next[i];
	}

	//	Compute gamma(D) and gamma'(D) in bulk
	(*Msm->sft->soften)(Msm->sft, Idx, D, F, DF);

	//	Compute energy and forces
	for (k = 0; k < Idx; k++)
	{
		//	U += q(i)*q(j)*(1/d - gamma(d/a)/a);
		(*U) += R[k].q*(1.0/(a*D[k]) - F[k]/a);

		//	f += q(i)*q(j)*(1/d^2 + gamma'(d/a)/a^2) * (r[j]-r[i])/d;
		Magnitude = R[k].q*(1.0/(a2*D[k]*D[k]) + DF[k]/a2);
		Direction[0] = R[k].x / (a*D[k]);
		Direction[1] = R[k].y / (a*D[k]);
		Direction[2] = R[k].z / (a*D[k]);

		f[I[k]][0] -= Magnitude*Direction[0];
		f[I[k]][1] -= Magnitude*Direction[1];
		f[I[k]][2] -= Magnitude*Direction[2];

		f[J[k]][0] += Magnitude*Direction[0];
		f[J[k]][1] += Magnitude*Direction[1];
		f[J[k]][2] += Magnitude*Direction[2];
	}

	//	Free dynamically allocated memory
	dynfree(D);
	dynfree(F);
	dynfree(DF);
	dynfree(R);
	dynfree(I);
	dynfree(J);
}

void msm_short_range_compute_neighbor(MSM* Msm, SIMULATION_DOMAIN* Domain, long* First, long* Next, long* ParticlesPerBin, long XBinCount, long YBinCount, long x, long y, long z, long l, long m, long n)
{
	long		i = -1;
	long		j = -1;
	long		k = 0;
	double		a = Msm->prm.a;
	double		a2 = a*a;
	PARTICLE*	r = Domain->Particles->r;
	double*		U = &Domain->Particles->U;
	double**	f = Domain->Particles->f;
	double		dx = 0.0;
	double		dy = 0.0;
	double		dz = 0.0;
	double		d = 0.0;
	long		Idx = 0;
	double*		D = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	PARTICLE*	R = NULL;
	long*		I = NULL;
	long*		J = NULL;
	double		Magnitude = 0.0;
	double		Direction[3];
	long		ParticlesInBin = 0;
	long		MaxInteractionCount = 0;

	//	What is the maximum number of interactions which we could compute?
	ParticlesInBin = ParticlesPerBin[IDX(x,y,z,XBinCount,YBinCount)];
	MaxInteractionCount = ParticlesInBin*ParticlesPerBin[IDX(l,m,n,XBinCount,YBinCount)];

	//	Dynamically allocate memory
	D = (double*) dynvec(MaxInteractionCount, sizeof(double));
	F = (double*) dynvec(MaxInteractionCount, sizeof(double));
	DF = (double*) dynvec(MaxInteractionCount, sizeof(double));
	R = (PARTICLE*) dynvec(MaxInteractionCount, sizeof(PARTICLE));
	I = (long*) dynvec(MaxInteractionCount, sizeof(long));
	J = (long*) dynvec(MaxInteractionCount, sizeof(long));

	//	Particle i in bin 1
	i = First[IDX(x,y,z,XBinCount,YBinCount)];
	while (i != -1)
	{
		j = First[IDX(l,m,n,XBinCount,YBinCount)];	//	Resets particle j to first in the list for new particle i

		//	Particle j in bin 2
		while (j != -1)
		{
			dx = r[j].x - r[i].x;
			dy = r[j].y - r[i].y;
			dz = r[j].z - r[i].z;
			d = sqrt(dx*dx + dy*dy + dz*dz);
			assert(d > 0.0);

			if (d < a)
			{
				//	Add d/a to vector D, NOTE: need a*D[x] to get d below
				D[Idx] = d/a;

				//	Add r[j]-r[i] direction vector and q[j]*q[i] to array R
				R[Idx].x = dx;
				R[Idx].y = dy;
				R[Idx].z = dz;
				R[Idx].q = r[j].q * r[i].q;

				//	Add i to vector I and add j to vector J
				I[Idx] = i;
				J[Idx] = j;

				Idx++;
			}
			j = Next[j];
		}
		i = Next[i];
	}

	//	Compute gamma(D) and gamma'(D) in bulk
	(*Msm->sft->soften)(Msm->sft, Idx, D, F, DF);

	//	Compute energy and forces
	for (k = 0; k < Idx; k++)
	{
		//	U += q(i)*q(j)*(1/d - gamma(d/a)/a);
		(*U) += R[k].q*(1.0/(a*D[k]) - F[k]/a);

		//	f += q(i)*q(j)*(1/d^2 + gamma'(d/a)/a^2) * (r[j]-r[i])/d;
		Magnitude = R[k].q*(1.0/(a2*D[k]*D[k]) + DF[k]/a2);
		Direction[0] = R[k].x / (a*D[k]);
		Direction[1] = R[k].y / (a*D[k]);
		Direction[2] = R[k].z / (a*D[k]);

		f[I[k]][0] -= Magnitude*Direction[0];
		f[I[k]][1] -= Magnitude*Direction[1];
		f[I[k]][2] -= Magnitude*Direction[2];

		f[J[k]][0] += Magnitude*Direction[0];
		f[J[k]][1] += Magnitude*Direction[1];
		f[J[k]][2] += Magnitude*Direction[2];
	}

	//	Free dynamically allocated memory
	dynfree(D);
	dynfree(F);
	dynfree(DF);
	dynfree(R);
	dynfree(I);
	dynfree(J);
}

//	INTERNAL TEST Methods
void msm_short_range_naive(MSM* Msm, SIMULATION_DOMAIN* Domain)
{
/*		-- In evaluate()
	double		U2 = 0.0;
	double**	f = Domain->Particles->f;
	double**	f2 = NULL;
	long		i = 0;
	double		maxferr = 0.0;

		--	In short_range if-clause
		//	Save off energy/forces for comparison
		U2 = Domain->Particles->U;
		f2 = (double**) dynarr_d(N,3);
		memcpy(f2[0], f[0], N*3*sizeof(double));

		//	Re-initialize output variables U and f
		Domain->Particles->U = 0.0;
		memset(f[0], 0, sizeof(double)*N*3);	//	FIXME <-- double check this

		msm_short_range_naive(Msm, Domain);

		//	Compare Energy and forces
		printf("\t\tDiff in Energy: %e\n", fabs(U2-Domain->Particles->U));

		for (i = 0; i < N; i++)
		{
			f2[i][0] = (f2[i][0]-f[i][0])*(f2[i][0]-f[i][0]);
			f2[i][1] = (f2[i][1]-f[i][1])*(f2[i][1]-f[i][1]);
			f2[i][2] = (f2[i][2]-f[i][2])*(f2[i][2]-f[i][2]);
			maxferr = MAX(maxferr, sqrt(f2[i][0] + f2[i][1] + f2[i][2]));
		}
		printf("\t\tDiff in force (max of N 2-norms): %e\n", maxferr);

		dynfree(f2[0]);
		dynfree(f2);
*/
	long		i = 0;
	long		j = 0;
	double		d = 0.0;
	double		dx = 0.0;
	double		dy = 0.0;
	double		dz = 0.0;
	double		a = Msm->prm.a;
	long		N = Domain->Particles->N;
	PARTICLE*	r = Domain->Particles->r;
	double*		U = &Domain->Particles->U;
	double**	f = Domain->Particles->f;
	double		dfx = 0.0;
	double		dfy = 0.0;
	double		dfz = 0.0;

	double		gamma = 0.0;
	double		dgamma = 0.0;
	double		d_a = 0.0;

	//	Perform the naive O(N^2) calculation
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			//	Compute Euclidean distance between two particles
			dx = r[j].x-r[i].x;
			dy = r[j].y-r[i].y;
			dz = r[j].z-r[i].z;
			d = sqrt(dx*dx + dy*dy + dz*dz);

			if (d < a)
			{
				d_a = d/a;
				//	Compute contribution to the energy
				(*Msm->sft->soften)(Msm->sft, 1, &d_a, &gamma, &dgamma);
				*U += (r[i].q*r[j].q*(1.0/d - gamma/a));

				//	Compute contribution to the forces
				dfx = (dx/d)*r[i].q*r[j].q*(1.0/(d*d) + dgamma/(a*a));
				dfy = (dy/d)*r[i].q*r[j].q*(1.0/(d*d) + dgamma/(a*a));
				dfz = (dz/d)*r[i].q*r[j].q*(1.0/(d*d) + dgamma/(a*a));

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
}

//	End of file
