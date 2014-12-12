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
void msm_bin_to_bin(MSM* Msm, SIMULATION_DOMAIN* Domain, long* Next, long Particle1, long Particle2, long MaxIterationCount)
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
	double*		D_A = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	PARTICLE*	R = NULL;
	long*		I = NULL;
	long*		J = NULL;
	double		Magnitude = 0.0;
	double		Direction[3];

	//	Dynamically allocate memory
	D = (double*) dynvec(MaxIterationCount, sizeof(double));
	D_A = (double*) dynvec(MaxIterationCount, sizeof(double));
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
				//	Add d to vector D
				D[Idx] = d;
				D_A[Idx] = d/a;

				//	Add q[i]*q[j] to vector Q	-> Save in PARTICLE array?
				//	Add r[j]-r[i] to array R	-> Save in PARTICLE array?
				R[Idx].x = dx;
				R[Idx].y = dy;
				R[Idx].z = dz;
				R[Idx].q = 1.0;//r[j].q * r[i].q;

				//	Add i to vector I
				//	add j to vector J
				I[Idx] = i;
				J[Idx] = j;
if (i == j) printf("\t\ti == j, OMG\n"); else printf("(i,j) = (%ld,%ld), d=%f, Idx=%ld\n", i,j,D[Idx], Idx);

				Idx++;
			}
			j = Next[j];
		}
		i = Next[i];
	}

//	if (Idx > MaxIterationCount) printf("Interactions=%ld, Max=%ld\n", Idx, MaxIterationCount);
printf("Interactions=%ld, Max=%ld\n", Idx, MaxIterationCount);

	//	Compute gamma(D) and gamma'(D) in bulk
	(*Msm->sft->soften)(Msm->sft, Idx, D_A, F, DF);

	//	Compute energy and forces
	for (k = 0; k < Idx; k++)
	{
//		printf("Idx=%ld\t%f\t%f\t%f\n", k, D[k], F[k], DF[k]);

		//	U += q(i)*q(j)*(1/d - gamma(d)/a);
		(*U) += R[k].q*(1.0/D[k] - F[k]/a);

		//	f += q(i)*q(j)*(1/d^2 + gamma'(d)/a^2) * (r[j]-r[i])/d;
		Magnitude = R[k].q*(1.0/(D[k]*D[k]) + DF[k]/a2);
		Direction[0] = R[k].x / D[k];
		Direction[1] = R[k].y / D[k];
		Direction[2] = R[k].z / D[k];

		f[I[k]][0] -= Magnitude*Direction[0];
		f[I[k]][1] -= Magnitude*Direction[1];
		f[I[k]][2] -= Magnitude*Direction[2];

		f[J[k]][0] += Magnitude*Direction[0];
		f[J[k]][1] += Magnitude*Direction[1];
		f[J[k]][2] += Magnitude*Direction[2];

	}

	//	Free dynamically allocated memory
	dynfree(D);
	dynfree(D_A);
	dynfree(F);
	dynfree(DF);
	dynfree(R);
	dynfree(I);
	dynfree(J);
}

void msm_short_range_compute_self(MSM* Msm, SIMULATION_DOMAIN* Domain, long* First, long* Next, long* ParticlesPerBin, long YBinCount, long ZBinCount, long i, long j, long k)
{
	long		ParticlesInBin = 0;
	long		MaxInteractionCount = 0;
	long		Particle1 = -1;
	long		Particle2 = -1;

	ParticlesInBin = ParticlesPerBin[IDX(i,j,k,YBinCount,ZBinCount)];
	MaxInteractionCount = ParticlesInBin*ParticlesInBin;
	if (MaxInteractionCount > 0)
	{
		Particle1 = First[IDX(i,j,k,YBinCount,ZBinCount)];
		//Particle2 = Next[Particle1];
		Particle2 = Particle1;
		msm_bin_to_bin(Msm, Domain, Next, Particle1, Particle2, MaxInteractionCount);
	}
}

void msm_short_range_compute_neighbor(MSM* Msm, SIMULATION_DOMAIN* Domain, long* First, long* Next, long* ParticlesPerBin, long YBinCount, long ZBinCount, long i, long j, long k, long l, long m, long n)
{
	long		ParticlesInBin = 0;
	long		MaxInteractionCount = 0;
	long		Particle1 = -1;
	long		Particle2 = -1;

	ParticlesInBin = ParticlesPerBin[IDX(i,j,k,YBinCount,ZBinCount)];
	MaxInteractionCount = ParticlesInBin*ParticlesPerBin[IDX(l,m,n,YBinCount,ZBinCount)];
	if (MaxInteractionCount > 0)
	{
		Particle1 = First[IDX(i,j,k,YBinCount,ZBinCount)];
		Particle2 = First[IDX(l,m,n,YBinCount,ZBinCount)];
		msm_bin_to_bin(Msm, Domain, Next, Particle1, Particle2, MaxInteractionCount);
	}
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
	long		Particle1 = 0;
	long		Particle2 = 0;
	long		MaxInteractionCount = 0;
	long*		First = NULL;
	long*		Next = NULL;
	long*		ParticlesPerBin = NULL;
	long		ParticlesInBin = 0;

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
		Next[n] = First[IDX(i,j,k,YBinCount,ZBinCount)];
		First[IDX(i,j,k,YBinCount,ZBinCount)] = n;
		ParticlesPerBin[IDX(i,j,k,YBinCount,ZBinCount)]++;
	}
printf("\t\tbin 0 has %ld particles in it\n", ParticlesPerBin[0]);
/*
	//	Loop over bins
	for (i = 0; i < XBinCount; i++)
	{
		for (j = 0; j < YBinCount; j++)
		{
			for (k = 0; k < ZBinCount; k++)
			{
				printf("%03ld <-- (%02ld,%02ld,%02ld), Cnt=%ld: ", IDX(i,j,k,YBinCount,ZBinCount), i,j,k, ParticlesPerBin[IDX(i,j,k,YBinCount,ZBinCount)]);
				if (First[IDX(i,j,k,YBinCount,ZBinCount)] != -1)
				{
//					printf("Bin: (%02ld,%02ld,%02ld):\t", i,j,k);
					n = First[IDX(i,j,k,YBinCount,ZBinCount)];
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

	//	i = 0 (no i-1)
	i = 0;
	{
		//	j = 0 (no j-1)
		j = 0;
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i+1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j,k-1);
				//	"next" bin is (i, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j+1,k-1);
				//	"next" bin is (i+1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j+1,k-1);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}

		for (j = 1; j < YBinCount-1; j++)
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k-1);
				//	"next" bin is (i+1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k-1);
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i+1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j,k-1);
				//	"next" bin is (i, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j+1,k-1);
				//	"next" bin is (i+1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j+1,k-1);
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}

		//	j = YBinCount (no j+1) (FIXME?: no j-1 if YBinCount == 1) NOTE: if YBinCount == 1, then j = 0 has already been computed above!
		if ((j = YBinCount-1) > 0)
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);//
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k-1);//
				//	"next" bin is (i+1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k-1);//
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i+1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j,k-1);
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);//
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}
	}

	for (i = 1; i < XBinCount-1; i++)
	{
		//	j = 0 (no j-1)
		j = 0;
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i-1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k-1);
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i+1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j,k-1);
				//	"next" bin is (i-1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j+1,k-1);
				//	"next" bin is (i, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j+1,k-1);
				//	"next" bin is (i+1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j+1,k-1);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}

		for (j = 1; j < YBinCount-1; j++)
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	ALL
				//	"next" bin is (i-1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k-1);
				//	"next" bin is (i, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k-1);
				//	"next" bin is (i+1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k-1);
				//	"next" bin is (i-1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k-1);
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i+1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j,k-1);
				//	"next" bin is (i-1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j+1,k-1);
				//	"next" bin is (i, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j+1,k-1);
				//	"next" bin is (i+1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j+1,k-1);
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}

		//	j = YBinCount (no j+1) (FIXME?: no j-1 if YBinCount == 1) NOTE: If YBinCount == 0, then j = 0 was already computed above!
		if ((j = YBinCount-1) > 0)
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);//
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);//
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);//
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i-1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k-1);//
				//	"next" bin is (i, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k-1);//
				//	"next" bin is (i+1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k-1);//
				//	"next" bin is (i-1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k-1);
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i+1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j,k-1);
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);//
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);//
				//	"next" bin is (i+1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i+1,j-1,k);//
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}
	}

	//	i = XBinCount (no i+1) (FIXME?: no i-1 if XBinCount == 1) NOTE: If XBinCount = 1, then i = 0 has already been computed above!
	if ((i = XBinCount-1) > 0)
	{
		//	j = 0 (no j-1)
		j = 0;
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i-1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k-1);//
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i-1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j+1,k-1);//
				//	"next" bin is (i, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j+1,k-1);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}

		for (j = 1; j < YBinCount-1; j++)
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);//
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i-1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k-1);//
				//	"next" bin is (i, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k-1);
				//	"next" bin is (i-1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k-1);//
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i-1, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j+1,k-1);//
				//	"next" bin is (i, j+1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j+1,k-1);
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);//
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}

		//	j = YBinCount (no j+1) (FIXME?: no j-1 if YBinCount == 1) NOTE: If YBinCount == 1 then j = 0 was already computed above!
		if ((j = YBinCount-1) > 0)
		{
			//	k = 0 (no k-1)
			k = 0;
			{
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);//
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}

			for (k = 1; k < ZBinCount; k++)
			{
				//	"next" bin is (i-1, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k-1);//
				//	"next" bin is (i, j-1, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k-1);
				//	"next" bin is (i-1, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k-1);//
				//	"next" bin is (i, j, k-1)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j,k-1);
				//	"next" bin is (i-1, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j-1,k);//
				//	"next" bin is (i, j-1, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i,j-1,k);
				//	"next" bin is (i-1, j, k)
				msm_short_range_compute_neighbor(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i,j,k, i-1,j,k);//
				//	"next" bin is (i, j, k)
				msm_short_range_compute_self(Msm, Domain, First, Next, ParticlesPerBin, YBinCount, ZBinCount, i, j, k);
			}
		}
	}

	//	Free dynamically allocated memory
	dynfree(First);
	dynfree(Next);
	dynfree(ParticlesPerBin);
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
