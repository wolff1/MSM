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
//	printf("Initializing MSM!\n");

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
	Msm->prm.L = 1;

	//	Initialize MSM options
	Msm->opt.ComputeExclusions = 0;
	Msm->opt.ComputeLongRange = 1;
	Msm->opt.ComputeShortRange = 0;
	Msm->opt.IsN = 0;
	Msm->opt.IsNLogN = 1;
	Msm->opt.GridType = 0;

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
//	printf("MSM Preprocessing for size <%f>!\n", DomainRadius);

	//	Compute Interpolant coefficients
	Msm->prm.D = DomainRadius;
	(*Msm->itp->compute_g2g)(Msm->itp, Msm->sft, &Msm->prm);
	(*Msm->itp->compute_tg2g)(Msm->itp, Msm->sft, &Msm->prm);
}

void msm_evaluate(void* Method, SIMULATION_DOMAIN* Domain)
{
	short		l = 0;	//	Lower case L, for level
	size_t		Size = 0;
	void		(*Init)(void*, SIMULATION_DOMAIN*, short) = NULL;
	GRID*		ChargeGrid = NULL;
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	assert(Domain != NULL);

	//	Initialize output variables U and f
	Domain->Particles->U = 0.0;
	memset(Domain->Particles->f[0], 0, sizeof(double)*Domain->Particles->N*3);

	if (Msm->opt.ComputeShortRange)
	{
		msm_short_range(Msm, Domain);
	}

	if (Msm->opt.ComputeLongRange)
	{
		//	SET UP CONTAINER FOR FINEST GRID
		switch (Msm->opt.GridType)
		{
		case 0:	//	RECTANGULAR_ROW_MAJOR_B_SPLINE Grid
		default:
			Size = sizeof(RECTANGULAR_ROW_MAJOR_B_SPLINE);
			Init = &rectangular_row_major_b_spline_initialize;
		}
		ChargeGrid = (GRID*) dynmem(Size);
		ChargeGrid->initialize = Init;

		//	COMPUTE LONG RANGE COMPONENT, O(N)
		if (Msm->opt.IsN)
		{
			msm_anterpolate(Msm, Domain, 0, ChargeGrid);

			for (l = 1; l < Msm->prm.L; l++)
			{	//	l is fine grid level
				msm_restrict(Msm);
				msm_direct(Msm);
			}

			msm_direct_top(Msm);

			for (l = Msm->prm.L-1, l > 0; l--)
			{	//	l is fine grid level
				msm_prolongate(Msm);
			}

			msm_interpolate(Msm, ChargeGrid);
		}

		//	COMPUTE LONG RANGE COMPONENT, O(N*log(N))
		if (Msm->opt.IsNLogN)
		{
			//	Intermediate Grid Level(s)
			for (l = 0; l < Msm->prm.L-1; l++)
			{
				msm_anterpolate(Msm, Domain, l, ChargeGrid);
				msm_direct(Msm);
				msm_interpolate(Msm, ChargeGrid);
			}
			//	Top Grid Level
			msm_anterpolate(Msm, Domain, Msm->prm.L-1, ChargeGrid);
			msm_direct_top(Msm);
			msm_interpolate(Msm, ChargeGrid);
		}

		if (Msm->opt.ComputeExclusions)
		{
			msm_exclude(Msm);
		}

		//	Free FINEST grid memory
		dynfree(ChargeGrid);
	}
}

void msm_uninitialize(void* Method)
{
	MSM*		Msm = (MSM*) Method;
	assert(Msm != NULL);
//	printf("Un-initializing MSM!\n");

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
//printf("\t\t%ldx%ldx%ld bins\n", XBinCount, YBinCount, ZBinCount);
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

void msm_anterpolate(MSM* Msm, SIMULATION_DOMAIN* Domain, short Level, GRID* Grid)
{
	//	INPUT:	particles
	//	OUTPUT:	finest grid
	PARTICLE*					r = Domain->Particles->r;
	PARTICLE*					Min = &Domain->MinimumCoordinates;
	short						p = Msm->prm.p;
	long						n = 0;
	double						h = 1.0;
	short						nu = 0;
	double						Dx = 0.0;
	double						Dy = 0.0;
	double						Dz = 0.0;
	double*						X = NULL;
	double*						FX = NULL;
	double*						DFX = NULL;
	double*						PhiX = NULL;
	double*						PhiY = NULL;
	double*						PhiZ = NULL;
	double*						dPhiX = NULL;//can remove after interpolation
	double*						dPhiY = NULL;//can remove after interpolation
	double*						dPhiZ = NULL;//can remove after interpolation
	double						ChargeZ = 0.0;
	double						ChargeYZ = 0.0;
	double						ChargeXYZ = 0.0;
	long						i = 0;
	long						j = 0;
	long						k = 0;
	long						x = 0;
	long						y = 0;
	long						z = 0;

	printf("\tMSM anterpolation!\n");

	//	Create Grid <Level>
	grid_initialize(Grid, Domain, Level, Msm->prm.h, Msm->prm.p);
	h = Grid->h;

	//	Initialize vector(s) used to hold interpolant input and output
	X = (double*) dynvec(3*p, sizeof(double));
	FX = (double*) dynvec(3*p, sizeof(double));
	DFX = (double*) dynvec(3*p, sizeof(double));

	PhiX = &FX[0];
	PhiY = &FX[p];
	PhiZ = &FX[2*p];

	dPhiX = &DFX[0];//can remove after interpolation
	dPhiY = &DFX[p];//can remove after interpolation
	dPhiZ = &DFX[2*p];//can remove after interpolation

	//	Loop through all particles	-> Spread particle charge onto grid in each dimension
	for (n = 0; n < Domain->Particles->N; n++)
	{
		//	Nearest grid point indices for particle
		Dx = (r[n].x-Min->x)/h;
		Dy = (r[n].y-Min->y)/h;
		Dz = (r[n].z-Min->z)/h;

		//	Transform (x,y,z) from nearest grid point to "lowest" grid point within support
		x = (long) floor(Dx) - (p>>1) + 1;
		y = (long) floor(Dy) - (p>>1) + 1;
		z = (long) floor(Dz) - (p>>1) + 1;

		Dx = Dx - x;
		Dy = Dy - y;
		Dz = Dz - z;

		//	Get distances to the nearest grid points
		for (nu = 0; nu < p; nu++)
		{
			X[nu] =		Dx - (double)nu;
			X[nu+p] =	Dy - (double)nu;
			X[nu+2*p] =	Dz - (double)nu;
		}

		//	Get interpolant values at the nearest grid points (in bulk)
		(*Msm->itp->evaluate)(Msm->itp, 3*p, X, FX, DFX);

//		for (nu = 0; nu < p; nu++)
//		{
//			printf("nu = %hd - (%+0f, %+0f, %+0f), (%+0f, %+0f, %+0f)\n", nu, X[nu], X[nu+p], X[nu+2*p], FX[nu], FX[nu+p], FX[nu+2*p]);
//		}

		//	Distribute point charge to grid
		for (k = 0; k < p; k++)
		{
			ChargeZ = PhiZ[k]*r[n].q;
			for (j = 0; j < p; j++)
			{
				ChargeYZ = PhiY[j]*ChargeZ;
				for (i = 0; i < p; i++)
				{
					ChargeXYZ = PhiX[i]*ChargeYZ;
					(*Grid->increment_grid_point_value)(Grid, x+i, y+j, z+k, ChargeXYZ);
				}
			}
		}
	}

	//	Display the grid as a sanity check
//	(*Grid->display)(Grid);

	//	Free dynamically allocated memory
	dynfree(DFX);
	dynfree(FX);
	dynfree(X);
}

void msm_restrict(MSM* Msm)
{
	//	INPUT:	fine grid
	//	OUTPUT:	coarse grid
//	printf("\tMSM restriction!\n");
/*
	all_ranges = *all fine grid points*
	for (m = 0; m < all_ranges.num; m++)
	{
		for (i = all_ranges[m].min; i <= all_ranges[m].max; i++)
		{
			phi_ranges = *corresponding coarse grid points* for fine grid point, i
			for (n = 0; n < phi_ranges.num; n++)
			{
				for (j = phi_ranges[n].min; j <= phi_ranges[n].max; j++)
				{
					coarse_grid[j] += Phi(g(i,j))*fine_grid[i];
				}
			}
		}
	}
*/
}

void msm_direct(MSM* Msm)
{
	//	INTPUT:	charge grid
	//	OUTPUT:	potential grid
//	printf("\tMSM direct computation!\n");
/*
	all_ranges = *all potential grid points*
	for (m = 0; m < all_ranges.num; m++)
	{
		for (i = all_ranges[m].min; i <= all_ranges[m].max; i++)
		{
			stencil_ranges = *corresponding charge grid points* for grid point, i
			for (n = 0; n < stencil_ranges.num; n++)
			{
				for (j = stencil_ranges[n].min; j <= stencil_ranges[n].max; j++)
				{
					potential_grid[i] += K(h(i,j))*charge_grid[j];
				}
			}
		}
	}
*/
}

void msm_direct_top(MSM* Msm)
{
	//	INTPUT:	charge grid
	//	OUTPUT:	potential grid
//	printf("\tMSM direct computation (top-level)!\n");
/*
	all_ranges = *all potential grid points*
	for (m = 0; m < all_ranges.num; m++)
	{
		for (i = all_ranges[m].min; i <= all_ranges[m].max; i++)
		{
			top_stencil_ranges = *corresponding charge grid points* for grid point, i (probably all charge grid points)
			for (n = 0; n < top_stencil_ranges.num; n++)
			{
				for (j = top_stencil_ranges[n].min; j <= top_stencil_ranges[n].max; j++)
				{
					potential_grid[i] += K_L(h(i,j))*charge_grid[j];
				}
			}
		}
	}
*/
}

void msm_prolongate(MSM* Msm)
{
	//	INPUT:	coarse grid
	//	OUTPUT:	fine grid
//	printf("\tMSM prolongation!\n");
/*
	all_ranges = *all fine grid points*
	for (m = 0; m < all_ranges.num; m++)
	{
		for (i = all_ranges[m].min; i <= all_ranges[m].max; i++)
		{
			phi_ranges = *corresponding coarse grid points* for fine grid point, i
			for (n = 0; n < phi_ranges.num; n++)
			{
				for (j = phi_ranges[n].min; j <= phi_ranges[n].max; j++)
				{
					fine_grid[i] += Phi(g(i,j))*coarse_grid[j];
				}
			}
		}
	}
*/
}

void msm_interpolate(MSM* Msm, GRID* ChargeGrid)
{
	//	INPUT:	finest potential gird, finest charge grid
	//	OUTPUT:	energy, forces
	printf("\tMSM interpolation!\n");
/*
	E = 0.5*q_0'*e_0;

	for (n = 0; n < N; n++)
	{
		i = f(r[n],h);
		f[n] += Phi(g(r[n],h,i)*potential_grid[i]
	}
*/
	grid_uninitialize(ChargeGrid);
}

void msm_exclude(MSM* Msm)
{
//	printf("\tMSM exclusions!\n");
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
