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
	Msm->prm.a = 2.0;//12.5
	Msm->prm.h = 1.0;//2.5;
	Msm->prm.alpha = Msm->prm.a / Msm->prm.h;
	Msm->prm.p = 4;
	Msm->prm.k = 4;
	Msm->prm.mu = 10;
	Msm->prm.D = 0.0;	//	Not known until preprocess/evaluate
	Msm->prm.L = 2;		//	# of grids

	//	Initialize MSM options
	Msm->opt.ComputeExclusions = 1;
	Msm->opt.ComputeLongRange = 1;
	Msm->opt.ComputeShortRange = 1;
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
	GRID**		ChargeGrid = NULL;
	GRID**		PotentialGrid = NULL;
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	assert(Domain != NULL);

	//	Initialize output variables U and f
	Domain->Particles->U = 0.0;
	memset(Domain->Particles->f[0], 0, 3*Domain->Particles->N*sizeof(double));

	if (Msm->opt.ComputeShortRange)
	{
		msm_short_range(Msm, Domain);
	}

	if (Msm->opt.ComputeLongRange)
	{
		//	SET UP GRID(S)
		switch (Msm->opt.GridType)
		{
		case 0:	//	RECTANGULAR_ROW_MAJOR_B_SPLINE Grid
		default:
			Size = sizeof(RECTANGULAR_ROW_MAJOR_B_SPLINE);
			Init = &rectangular_row_major_b_spline_initialize;
		}
		//	Create and set initialization routine for array of anterpolation/restriction hierarchy CHARGE grids
		ChargeGrid = (GRID**) dynmem(sizeof(GRID*)*Msm->prm.L);
		PotentialGrid = (GRID**) dynmem(sizeof(GRID*)*Msm->prm.L);	//	FIXME: !!!This can be implemented with only 2 grid containers!!!
		for (l = 0; l < Msm->prm.L; l++)
		{
			ChargeGrid[l] = (GRID*) dynmem(Size);
			PotentialGrid[l] = (GRID*) dynmem(Size);
			((GRID*)ChargeGrid[l])->initialize = Init;
			((GRID*)PotentialGrid[l])->initialize = Init;
		}

		//	COMPUTE LONG RANGE COMPONENT, O(N)
		if (Msm->opt.IsN)
		{
			msm_anterpolate(Msm, Domain, 0, /*OUT*/ChargeGrid[0]);

			for (l = 1; l < Msm->prm.L; l++)	//	l is fine grid level
			{
				msm_restrict(Msm, /*IN*/ChargeGrid[l-1], /*OUT*/ChargeGrid[l]);
			}

			msm_direct_top(Msm, /*IN*/ChargeGrid[Msm->prm.L-1], /*OUT*/PotentialGrid[Msm->prm.L-1]);

			for (l = Msm->prm.L-1; l > 0; l--)	//	l is fine grid level
			{
				msm_direct(Msm, /*IN*/ChargeGrid[l-1], /*OUT*/PotentialGrid[l-1]);
				msm_prolongate(Msm, /*OUT*/PotentialGrid[l-1], /*IN*/PotentialGrid[l]);
			}

			msm_interpolate(Msm, Domain, /*IN*/ChargeGrid[0], /*IN*/PotentialGrid[0]);
		}

		//	COMPUTE LONG RANGE COMPONENT, O(N*log(N))
		if (Msm->opt.IsNLogN)
		{
//FIXME: Save off U and f in case O(N) method was also performed. U and f need to be re-initialized before proceeding here!
			//	Intermediate Grid Level(s)
			for (l = 0; l < Msm->prm.L-1; l++)
			{
				msm_anterpolate(Msm, Domain, l, /*OUT*/ChargeGrid[l]);
				msm_direct(Msm, /*IN*/ChargeGrid[l], /*OUT*/PotentialGrid[l]);
				msm_interpolate(Msm, Domain, /*IN*/ChargeGrid[l], /*IN*/PotentialGrid[l]);
			}
			//	Top Grid Level
			msm_anterpolate(Msm, Domain, Msm->prm.L-1, /*OUT*/ChargeGrid[Msm->prm.L-1]);
			msm_direct_top(Msm, /*IN*/ChargeGrid[Msm->prm.L-1], /*OUT*/PotentialGrid[Msm->prm.L-1]);
			msm_interpolate(Msm, Domain, /*IN*/ChargeGrid[Msm->prm.L-1], /*IN*/PotentialGrid[Msm->prm.L-1]);
		}

		if (Msm->opt.ComputeExclusions)
		{
			msm_exclude(Msm, Domain);
		}

		//	Free individual grid pointers
		for (l = 0; l < Msm->prm.L; l++)
		{
			dynfree(ChargeGrid[l]);
			dynfree(PotentialGrid[l]);
		}

		//	Free array of grid pointers
		dynfree(ChargeGrid);
		dynfree(PotentialGrid);
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
	long		MyIdx = 0;

	printf("\tMSM short range computation!\n");

	//	Set up number of bins in each dimension (NOTE: MATLAB code has floor() + 1)
	XBinCount = (long) ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / a) + 1;	// NOTE: +1 b/c if r[n] is at max, then it is in bin # XBinCount, NOT XBinCount-1
	YBinCount = (long) ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / a) + 1;	//			99% of the time, these +1 will NOT be necessary
	ZBinCount = (long) ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / a) + 1;
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
		MyIdx = IDX(i,j,k,XBinCount,YBinCount);
//printf("(%ld,%ld,%ld) = %ld\n", i,j,k, MyIdx);
		Next[n] = First[MyIdx];
		First[MyIdx] = n;
		ParticlesPerBin[MyIdx]++;
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
	long						i = 0;
	long						j = 0;
	long						k = 0;
	long						Di = 0;
	long						Dj = 0;
	long						Dk = 0;
	long						Idx = 0;
	long						ThisIdx = 0;
	double						ThisValue = 0.0;

	printf("\tMSM anterpolation!\n");

	//	Create Grid <Level> --> Freed in interpolation
	grid_initialize(Grid, Domain, Level, Msm->prm.h, Msm->prm.p);
printf("after grid initialize in anterpolate!\n");
	h = Grid->h;

	//	Initialize vector(s) used to hold interpolant input and output
	X = (double*) dynvec(3*p, sizeof(double));
	FX = (double*) dynvec(3*p, sizeof(double));
	DFX = (double*) dynvec(3*p, sizeof(double));

	PhiX = &FX[0];
	PhiY = &FX[p];
	PhiZ = &FX[2*p];

	//	Loop through all particles	-> Spread particle charge onto grid in each dimension
	for (n = 0; n < Domain->Particles->N; n++)
	{
		//	Nearest grid point indices for particle
		Dx = (r[n].x-Min->x)/h;
		Dy = (r[n].y-Min->y)/h;
		Dz = (r[n].z-Min->z)/h;

		//	Transform (x,y,z) from nearest grid point to "lowest" grid point within support
		i = (long) floor(Dx) - (p>>1) + 1;
		j = (long) floor(Dy) - (p>>1) + 1;
		k = (long) floor(Dz) - (p>>1) + 1;

		Dx = Dx - i;
		Dy = Dy - j;
		Dz = Dz - k;

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

		////	Distribute point charge to grid
		//for (Dk = 0; Dk < p; Dk++)
		//{
		//	ChargeZ = PhiZ[Dk]*r[n].q;
		//	for (Dj = 0; Dj < p; Dj++)
		//	{
		//		ChargeYZ = PhiY[Dj]*ChargeZ;
		//		for (Di = 0; Di < p; Di++)
		//		{
		//			ChargeXYZ = PhiX[Di]*ChargeYZ;
		//			(*Grid->increment_grid_point_value)(Grid, i+Di, j+Dj, k+Dk, ChargeXYZ);
		//		}
		//	}
		//}

		//	Distribute point charge to grid (IS THIS "BETTER" THAN ABOVE B/C LOOP IS "UNROLLED"?)
		for (Idx = 0; Idx < p*p*p; Idx++)
		{
			Dk = Idx / (p*p);
			Dj = (Idx - Dk*p*p) / p;
			Di = Idx - Dk*p*p - Dj*p;
//printf("(%ld,%ld,%ld) -> %ld\n", Di,Dj,Dk,Idx);
			ThisIdx = (*Grid->ijk2idx)(Grid, i+Di, j+Dj, k+Dk);
			ThisValue = PhiX[Di]*PhiY[Dj]*PhiZ[Dk]*r[n].q;
			(*Grid->increment_grid_point_value)(Grid, ThisIdx, ThisValue);
		}
	}

	//	Display the grid as a sanity check
//	(*Grid->display)(Grid);

	//	Free dynamically allocated memory
	dynfree(DFX);
	dynfree(FX);
	dynfree(X);
}

void msm_restrict(MSM* Msm, GRID* FineGrid, GRID* CoarseGrid)
{
	//	INPUT:	fine grid
	//	OUTPUT:	coarse grid

	//	Create Coarse Grid for <Level> --> Freed either in direct or direct_top
//	grid_initialize(Grid, Domain, Level, Msm->prm.h, Msm->prm.p);
//	grid_create_coarse_structure(?);

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

void msm_direct(MSM* Msm, GRID* ChargeGrid, GRID* PotentialGrid)
{
	//	INTPUT:	charge grid
	//	OUTPUT:	potential grid

	GRID_RANGE			Range;
	long				MaxSlices = 0;

	printf("\tMSM direct computation! Stencil Radius: %ld\n", Msm->itp->g2g->Size);

	//	Ranges is an array of length NumSlices
	MaxSlices = (2*(Msm->itp->g2g->Size)+1)*(2*(Msm->itp->g2g->Size)+1);
	Range.Ranges = (GRID_RANGE_MIN_MAX*) dynvec(MaxSlices,sizeof(GRID_RANGE_MIN_MAX));

	rectangular_row_major_b_spline_get_grid_points_stencil(ChargeGrid, rectangular_row_major_b_spline_ijk2idx(ChargeGrid,0,0,0), Msm->itp->g2g, &Range);

//	for (MaxSlices = 0; MaxSlices < Range.NumSlices; MaxSlices++)
//	{
//		printf("Slice %03ld: (%04ld, %04ld)\n", MaxSlices, Range.Ranges[MaxSlices].Min, Range.Ranges[MaxSlices].Max);
//	}

	//	Create Potential Grid for <Level> --> Freed either in interpolate or prolongation
	(*ChargeGrid->create_copy_grid_structure)(PotentialGrid, ChargeGrid);
//printf("After grid structure copy!\n");

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
	//	Free dynamically allocated memory
	dynfree(Range.Ranges);

	//	Free Charge Grid if not finest charge grid (b/c its needed in interpolation)
	if ((Msm->opt.IsN) && (ChargeGrid->Level > 0))
	{
		grid_uninitialize(ChargeGrid);
	}
}

void msm_direct_top(MSM* Msm, GRID* ChargeGrid, GRID* PotentialGrid)
{
	//	INTPUT:	charge grid
	//	OUTPUT:	potential grid

	GRID_RANGE		Outer;
	GRID_RANGE		Inner;
	long			MaxSlices = 0;
	long			m = 0;
	long			n = 0;
	long			i = 0;
	long			j = 0;
	double			GridValue = 0.0;
	long			Idx = 0;
	long			i1 = 0;
	long			j1 = 0;
	long			k1 = 0;
	long			i2 = 0;
	long			j2 = 0;
	long			k2 = 0;
	long			di = 0;
	long			dj = 0;
	long			dk = 0;
	long			MaxIdx = 0;
	long			MinIdx = 220;
	double			a_l = 1.0;

	printf("\tMSM direct computation (top-level)!\n");

	for (i = 0; i < ChargeGrid->Level; i++)
	{
		//	NOTE: Actual a_l = (2^Level)*a, this is 1/a_l
		a_l *= 0.5;
	}
	a_l /= Msm->prm.a;

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
	//	Create Potential Grid for <Level> --> Freed either in interpolate or prolongation
	(*ChargeGrid->create_copy_grid_structure)(PotentialGrid, ChargeGrid);

	//	Ranges is an array of length NumSlices
	MaxSlices = (*ChargeGrid->get_grid_points_all_max_slices)(ChargeGrid);
	Outer.Ranges = (GRID_RANGE_MIN_MAX*) dynvec(MaxSlices,sizeof(GRID_RANGE_MIN_MAX));
	Inner.Ranges = (GRID_RANGE_MIN_MAX*) dynvec(MaxSlices,sizeof(GRID_RANGE_MIN_MAX));//FIXME - use value from STENCIL instead of MaxSlices?

	//	Loop over *ALL* grid points to fill in potential grid from charge grid
	(*ChargeGrid->get_grid_points_all)(ChargeGrid, &Outer);
//	printf("Outer grid has <%ld> slices:\n", Outer.NumSlices);
	for (m = 0; m < Outer.NumSlices; m++)
	{
//		printf("Slice <%04ld> - [%04ld, %04ld]\n", m, Outer.Ranges[m].Min, Outer.Ranges[m].Max);
		for (i = Outer.Ranges[m].Min; i <= Outer.Ranges[m].Max; i++)
		{
			//	Loop over *TOP STENCIL* grid points to interact with those from the outer loop
			(*ChargeGrid->get_grid_points_stencil_top)(ChargeGrid, i, &Inner);
//			printf("\tInner grid has <%ld> slices:\n", Inner.NumSlices);
			for (n = 0; n < Inner.NumSlices; n++)
			{
//				printf("\tSlice <%04ld> - [%04ld, %04ld]\n", n, Inner.Ranges[n].Min, Inner.Ranges[n].Max);
				for (j = Inner.Ranges[n].Min; j <= Inner.Ranges[n].Max; j++)
				{
//					PotentialGrid[i] += K[h(i,j)]*ChargeGrid[j];
/*
					i -> (x1,y1,z1);
					j -> (x2,y2,z2);
*/
					(*ChargeGrid->idx2ijk)(ChargeGrid, i, &i1, &j1, &k1);
					(*ChargeGrid->idx2ijk)(ChargeGrid, j, &i2, &j2, &k2);
//FIXME - The following is cringe-worthy:
					di = abs(i2-i1);
					dj = abs(j2-j1);
					dk = abs(k2-k1);
					if (di <= dj)
					{
						if (dj <= dk)
						{
							//	i <= j <= k
							Idx = STENCIL_MAP_X(di) + STENCIL_MAP_Y(dj) + STENCIL_MAP_Z(dk);
						}
						else
						{
							if (di <= dk)
							{
								//	i <= k <= j
								Idx = STENCIL_MAP_X(di) + STENCIL_MAP_Y(dk) + STENCIL_MAP_Z(dj);
							}
							else
							{
								//	k <= i <= j
								Idx = STENCIL_MAP_X(dk) + STENCIL_MAP_Y(di) + STENCIL_MAP_Z(dj);
							}
						}
					}
					else
					{
						if (dk <= dj)
						{
							//	k <= j <= i
							Idx = STENCIL_MAP_X(dk) + STENCIL_MAP_Y(dj) + STENCIL_MAP_Z(di);
						}
						else
						{
							if (di <= dk)
							{
								//	j <= i <= k
								Idx = STENCIL_MAP_X(dj) + STENCIL_MAP_Y(di) + STENCIL_MAP_Z(dk);
							}
							else
							{
								//	j <= k <= i
								Idx = STENCIL_MAP_X(dj) + STENCIL_MAP_Y(dk) + STENCIL_MAP_Z(di);
							}
						}
					}
					if (Idx > MaxIdx)
						MaxIdx = Idx;
					if (Idx < MinIdx)
						MinIdx = Idx;
//printf("i=%+04ld (%+04ld,%+04ld,%+04ld), j=%+04ld (%+04ld,%+04ld,%+04ld), Idx=%+04ld, K=%f\n", i,i1,j1,k1, j,i2,j2,k2, Idx, Msm->itp->tg2g->Data[Idx]);
//	NOTE: stencil access should be more encapsulated!
					GridValue = (Msm->itp->tg2g->Data[Idx])*(*ChargeGrid->get_grid_point_value)(ChargeGrid, j);
					(*PotentialGrid->increment_grid_point_value)(PotentialGrid, i, GridValue*a_l);
				}
			}
		}
	}
	printf("Min/Max = %ld/%ld\n", MinIdx, MaxIdx);
//	(*PotentialGrid->display)(PotentialGrid);

	//	Free dynamically allocated memory
	dynfree(Inner.Ranges);
	dynfree(Outer.Ranges);

	//	Free Charge Grid if not finest charge grid (b/c its needed in interpolation)
	if ((Msm->opt.IsN) && (ChargeGrid->Level > 0))
	{
		grid_uninitialize(ChargeGrid);
	}
}

void msm_prolongate(MSM* Msm, GRID* FineGrid, GRID* CoarseGrid)
{
	//	INPUT:	coarse grid
	//	OUTPUT:	fine grid
//	printf("\tMSM prolongation!\n");

	//	Create Fine Grid for <Level> --> Freed either in interpolate or prolongation
//	grid_initialize(Grid, Domain, Level, Msm->prm.h, Msm->prm.p);
//	grid_create_fine_structure(?);

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
	//	Free Coarse Grid
	grid_uninitialize(CoarseGrid);
}

void msm_interpolate(MSM* Msm, SIMULATION_DOMAIN* Domain, GRID* ChargeGrid, GRID* PotentialGrid)
{
	//	INPUT:	finest potential gird, finest charge grid
	//	OUTPUT:	energy, forces
	GRID_RANGE					Outer;
	long						MaxSlices = 0;
	long						m = 0;
	double						Energy = 0.0;

	PARTICLE*					r = Domain->Particles->r;
	PARTICLE*					Min = &Domain->MinimumCoordinates;
	short						p = Msm->prm.p;
	long						n = 0;
	double						h = ChargeGrid->h;
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
	double*						dPhiX = NULL;
	double*						dPhiY = NULL;
	double*						dPhiZ = NULL;
	long						i = 0;
	long						j = 0;
	long						k = 0;
	long						Di = 0;
	long						Dj = 0;
	long						Dk = 0;
	long						Idx = 0;
	long						GridIdx = 0;
	double						GridValue = 0.0;
	double						Fx = 0.0;
	double						Fy = 0.0;
	double						Fz = 0.0;

	printf("\tMSM interpolation! Level<%hd>\n", ChargeGrid->Level);
/*
	E = 0.5*q_0'*e_0;

	for (n = 0; n < N; n++)
	{
		i = f(r[n],h);
		f[n] += Phi(g(r[n],h,i)*potential_grid[i]
	}
*/

//	*** FIRST, CALCULATE THE ELECTROSTATIC ENERGY ***
//FIXME: Consider doing this calculation in prolongation. Instead of storing e_0, do the dot product with q_0 as e_0 is formed.
	//	Ranges is an array of length NumSlices
	MaxSlices = (*ChargeGrid->get_grid_points_all_max_slices)(ChargeGrid);
	Outer.Ranges = (GRID_RANGE_MIN_MAX*) dynvec(MaxSlices,sizeof(GRID_RANGE_MIN_MAX));
	//	Loop over *ALL* grid points
	(*ChargeGrid->get_grid_points_all)(ChargeGrid, &Outer);
	for (m = 0; m < Outer.NumSlices; m++)
	{
		for (i = Outer.Ranges[m].Min; i <= Outer.Ranges[m].Max; i++)
		{
			Energy += (*ChargeGrid->get_grid_point_value)(ChargeGrid, i)*(*PotentialGrid->get_grid_point_value)(PotentialGrid, i);
		}
	}
	Domain->Particles->U += 0.5*Energy;

	//	Free dynamically allocated memory
	grid_uninitialize(ChargeGrid);
	dynfree(Outer.Ranges);

//	*** SECOND, CALCULATE THE CORRESPONDING FORCES ***
	//	Initialize vector(s) used to hold interpolant input and output
	X = (double*) dynvec(3*p, sizeof(double));
	FX = (double*) dynvec(3*p, sizeof(double));
	DFX = (double*) dynvec(3*p, sizeof(double));

	PhiX = &FX[0];
	PhiY = &FX[p];
	PhiZ = &FX[2*p];

	dPhiX = &DFX[0];
	dPhiY = &DFX[p];
	dPhiZ = &DFX[2*p];

	//	Loop through all particles	-> Use grid potential to interpolate inter-particle forces
	for (n = 0; n < Domain->Particles->N; n++)
	{
		//	Nearest grid point indices for particle
		Dx = (r[n].x-Min->x)/h;
		Dy = (r[n].y-Min->y)/h;
		Dz = (r[n].z-Min->z)/h;

		//	Transform (x,y,z) from nearest grid point to "lowest" grid point within support
		i = (long) floor(Dx) - (p>>1) + 1;
		j = (long) floor(Dy) - (p>>1) + 1;
		k = (long) floor(Dz) - (p>>1) + 1;

		Dx = Dx - i;
		Dy = Dy - j;
		Dz = Dz - k;

		//	Get distances to the nearest grid points
		for (nu = 0; nu < p; nu++)
		{
			X[nu] =		Dx - (double)nu;
			X[nu+p] =	Dy - (double)nu;
			X[nu+2*p] =	Dz - (double)nu;
		}

		//	Get interpolant values at the nearest grid points (in bulk)
		(*Msm->itp->evaluate)(Msm->itp, 3*p, X, FX, DFX);

		//	Gather contributions to particle n forces
		Fx = 0.0;
		Fy = 0.0;
		Fz = 0.0;

		//	Gather contributions to particle n forces (IS THIS "BETTER" THAN ABOVE B/C LOOP IS "UNROLLED"?)
		for (Idx = 0; Idx < p*p*p; Idx++)
		{
			Dk = Idx / (p*p);
			Dj = (Idx - Dk*p*p) / p;
			Di = Idx - Dk*p*p - Dj*p;
//printf("(%ld,%ld,%ld) -> %ld\n", Di,Dj,Dk,Idx);
			GridIdx = (*PotentialGrid->ijk2idx)(PotentialGrid, i+Di, j+Dj, k+Dk);
			GridValue = (*PotentialGrid->get_grid_point_value)(PotentialGrid, GridIdx);
			Fx += dPhiX[Di]* PhiY[Dj]* PhiZ[Dk]*GridValue;
			Fy +=  PhiX[Di]*dPhiY[Dj]* PhiZ[Dk]*GridValue;
			Fz +=  PhiX[Di]* PhiY[Dj]*dPhiZ[Dk]*GridValue;
		}
		Domain->Particles->f[n][0] -= Fx*r[n].q/h;
		Domain->Particles->f[n][1] -= Fy*r[n].q/h;
		Domain->Particles->f[n][2] -= Fz*r[n].q/h;
	}

	//	Free more dynamically allocated memory
	grid_uninitialize(PotentialGrid);
	dynfree(DFX);
	dynfree(FX);
	dynfree(X);
}

void msm_exclude(MSM* Msm, SIMULATION_DOMAIN* Domain)
{
//	printf("\tMSM exclusions!\n");
	double		X = 0.0;
	double		FX = 0.0;
	double		DFX = 0.0;
	double		SelfEnergy = 0.0;
	long		i = 0;

	//	Self Energy = (1/2) (q^T,q) gamma(0.0)
	(*Msm->sft->soften)(Msm->sft, 1, &X, &FX, &DFX);

	for (i = 0; i < Domain->Particles->N; i++)
	{
		SelfEnergy += Domain->Particles->r[i].q*Domain->Particles->r[i].q;
	}

	Domain->Particles->U -= 0.5*FX*SelfEnergy/Msm->prm.a;
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
