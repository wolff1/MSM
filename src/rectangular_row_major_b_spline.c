//-------|---------|---------|---------|---------|---------|---------|---------|
/*
rectangular_row_major_b_spline.c -
*/

#include "rectangular_row_major_b_spline.h"

//	EXTERNAL Methods
void	rectangular_row_major_b_spline_initialize(void* Grid, SIMULATION_DOMAIN* Domain, short Expansion /*, alpha, stencil_shape*/)
{
	RECTANGULAR_ROW_MAJOR_B_SPLINE*		MyGrid = (RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid;

	assert(MyGrid != NULL);

//	printf("Min(%4.2f,%4.2f,%4.2f) Max(%4.2f,%4.2f,%4.2f) Len(%4.2f,%4.2f,%4.2f)\n",
//		Domain->MinimumCoordinates.x,Domain->MinimumCoordinates.y,Domain->MinimumCoordinates.z,
//		Domain->MaximumCoordinates.x,Domain->MaximumCoordinates.y,Domain->MaximumCoordinates.z,
//		Domain->MaximumCoordinates.x-Domain->MinimumCoordinates.x,Domain->MaximumCoordinates.y-Domain->MinimumCoordinates.y,Domain->MaximumCoordinates.z-Domain->MinimumCoordinates.z);

	//	Initialize COMMON members

	//	Initialize COMMON function pointers
	MyGrid->cmn.copy							= &rectangular_row_major_b_spline_copy;
	MyGrid->cmn.display							= &rectangular_row_major_b_spline_display;
	MyGrid->cmn.xyz2idx							= &rectangular_row_major_b_spline_xyz2idx;
	MyGrid->cmn.ijk2idx							= &rectangular_row_major_b_spline_ijk2idx;
	MyGrid->cmn.idx2ijk							= &rectangular_row_major_b_spline_idx2ijk;
	MyGrid->cmn.get_grid_points_all_max_slices	= &rectangular_row_major_b_spline_get_grid_points_all_max_slices;
	MyGrid->cmn.get_grid_points_all				= &rectangular_row_major_b_spline_get_grid_points_all;
	MyGrid->cmn.get_grid_points_coarse			= &rectangular_row_major_b_spline_get_grid_points_coarse;
	MyGrid->cmn.get_grid_points_stencil			= &rectangular_row_major_b_spline_get_grid_points_stencil;
	MyGrid->cmn.get_grid_points_stencil_top		= &rectangular_row_major_b_spline_get_grid_points_stencil_top;
	MyGrid->cmn.get_grid_point_value			= &rectangular_row_major_b_spline_get_grid_point_value;
	MyGrid->cmn.increment_grid_point_value		= &rectangular_row_major_b_spline_increment_grid_point_value;
	MyGrid->cmn.create_copy_grid_structure		= &rectangular_row_major_b_spline_create_copy_grid_structure;
	MyGrid->cmn.create_finer_grid				= &rectangular_row_major_b_spline_create_finer_grid;
	MyGrid->cmn.create_coarser_grid				= &rectangular_row_major_b_spline_create_coarser_grid;
	MyGrid->cmn.uninitialize					= &rectangular_row_major_b_spline_uninitialize;

	//	Initialize RECTANGULAR_ROW_MAJOR_B_SPLINE members
	MyGrid->Expansion = Expansion;
	MyGrid->Nx = (long) 2*ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / (2*MyGrid->cmn.h)) + Expansion + 1;
	MyGrid->Ny = (long) 2*ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / (2*MyGrid->cmn.h)) + Expansion + 1;
	MyGrid->Nz = (long) 2*ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / (2*MyGrid->cmn.h)) + Expansion + 1;
//printf("Grid Level <%ld> (%ld, %ld, %ld)\n", MyGrid->cmn.Level, MyGrid->Nx, MyGrid->Ny, MyGrid->Nz);

	MyGrid->Data = (double*) dynvec(MyGrid->Nx*MyGrid->Ny*MyGrid->Nz, sizeof(double));

	MyGrid->NumIntergridCoefficients = Expansion/2 + 1;
}

void	rectangular_row_major_b_spline_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	Create copy of structure
	rectangular_row_major_b_spline_create_copy_grid_structure(Dst, Src);

	//	Copy data from Src to Dst
	//		-> FIXME
}

void	rectangular_row_major_b_spline_display(void* Grid)
{
	long								i = 0;
	long								j = 0;
	long								k = 0;
	long								Idx = 0;
	long								q = 0;
	long								r = 0;
	long								s = 0;
	RECTANGULAR_ROW_MAJOR_B_SPLINE*		MyGrid = (RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid;
	short								Exp = MyGrid->Expansion >> 1;

	printf("\nRectangular, Row-major, B-spline Grid %hd: (%ld,%ld,%ld) -> %ld:\n\n",
			MyGrid->cmn.Level, MyGrid->Nx, MyGrid->Ny, MyGrid->Nz, MyGrid->Nx*MyGrid->Ny*MyGrid->Nz);

	for (k = 0; k < MyGrid->Nz; k++)
	{
		for (i = 0; i < MyGrid->Nx; i++)
		{
			for (j = 0; j < MyGrid->Ny; j++)
			{
				printf("%+06.3f ", MyGrid->Data[IDX(i,j,k,MyGrid->Nx,MyGrid->Ny)]);
			}
			printf("\n");
		}
		printf("\n");
	}

/*
	for (k = 0; k < MyGrid->Nz; k++)
	{
		for (j = 0; j < MyGrid->Ny; j++)
		{
			for (i = 0; i < MyGrid->Nx; i++)
			{
				Idx = (*MyGrid->cmn.ijk2idx)(Grid, i-Exp,j-Exp,k-Exp);
				(*MyGrid->cmn.idx2ijk)(Grid, Idx, &q, &r, &s);
				printf("(%03ld,%03ld,%03ld) -> %ld -> (%03ld,%03ld,%03ld)\n", i-Exp,j-Exp,k-Exp, Idx, q,r,s);
			}
			printf("\n");
		}
		printf("\n");
	}
*/
}

long	rectangular_row_major_b_spline_xyz2idx(void* Grid, SIMULATION_DOMAIN* Domain, double x, double y , double z)
{
	short				Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Expansion;
	double				h = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->cmn.h;
	PARTICLE*			Min = &Domain->MinimumCoordinates;

	return rectangular_row_major_b_spline_ijk2idx(Grid,
												  (long)floor((x - Min->x)/h),
												  (long)floor((y - Min->y)/h),
												  (long)floor((z - Min->z)/h));
}

long	rectangular_row_major_b_spline_ijk2idx(void* Grid, long i, long j, long k)
{
	short				Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Expansion >> 1;
	long				Nx =  ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Nx;
	long				Ny =  ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Ny;

	return IDX(i+Exp, j+Exp, k+Exp, Nx, Ny);
}

void	rectangular_row_major_b_spline_idx2ijk(void* Grid, long Idx, long* i, long* j, long* k)
{
	short				Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Expansion >> 1;
	long				Nx = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Nx;
	long				Ny = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Ny;
	long				Nz = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Nz;

	//	Convert Idx into components (i,j,k)
	*k = Idx / (Nx*Ny) - Exp;
	*j = (Idx - (*k+Exp)*Nx*Ny) / Nx - Exp;
	*i = Idx - (*k+Exp)*Nx*Ny - (*j+Exp)*Nx - Exp;
}

long	rectangular_row_major_b_spline_get_grid_points_all_max_slices(void* Grid)
{
	RECTANGULAR_ROW_MAJOR_B_SPLINE*		MyGrid = (RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid;
	return MyGrid->Ny*MyGrid->Nz;
}

void	rectangular_row_major_b_spline_get_grid_points_all(void* Grid, GRID_RANGE* Range)
{
	//	Return grid point indices for entire grid, from low to high
	//		indices are memory array indices, not (i,j,k)
	//		-> Should it be this way?
	long								j = 0;
	long								k = 0;
	RECTANGULAR_ROW_MAJOR_B_SPLINE*		MyGrid = (RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid;

	Range->NumSlices = 0;
	for (k = 0; k < MyGrid->Nz; k++)
	{
		for (j = 0; j < MyGrid->Ny; j++)
		{
			Range->Ranges[Range->NumSlices].Min = k*MyGrid->Nx*MyGrid->Ny + j*MyGrid->Nx;
			Range->Ranges[Range->NumSlices].Max = Range->Ranges[Range->NumSlices].Min + MyGrid->Nx-1;
			Range->NumSlices++;
		}
	}
}

void	rectangular_row_major_b_spline_get_grid_points_coarse(void* FineGrid, long FineGridIndex, void* CoarseGrid, GRID_RANGE* Range)
{
	long				i = 0;
	long				j = 0;
	long				k = 0;
	long				y = 0;
	long				z = 0;
	long				iMin = 0;
	long				jMin = 0;
	long				kMin = 0;
	long				iMax = 0;
	long				jMax = 0;
	long				kMax = 0;

	//	FineGridIndex -> (i,j,k)
	rectangular_row_major_b_spline_idx2ijk(FineGrid, FineGridIndex, &i, &j, &k);
	iMin = (long)ceil(i/2.0)-1;
	jMin = (long)ceil(j/2.0)-1;
	kMin = (long)ceil(k/2.0)-1;
	iMax = (long)floor(i/2.0)+1;
	jMax = (long)floor(j/2.0)+1;
	kMax = (long)floor(k/2.0)+1;
//printf("Fine(%ld,%ld,%ld), Coarse(%ld,%ld,%ld) -> (%ld,%ld,%ld)\n", i,j,k, iMin,jMin,kMin, iMax,jMax,kMax);

	Range->NumSlices = 0;
	for (z = kMin; z <= kMax; z++)
	{
		for (y = jMin; y <= jMax; y++)
		{
			Range->Ranges[Range->NumSlices].Min = rectangular_row_major_b_spline_ijk2idx(CoarseGrid, iMin,y,z);
			Range->Ranges[Range->NumSlices].Max = rectangular_row_major_b_spline_ijk2idx(CoarseGrid, iMax,y,z);
			Range->NumSlices++;
		}
	}
//printf("Num slices: %ld\n", Range->NumSlices);
}

void	rectangular_row_major_b_spline_get_grid_points_stencil(void* Grid, long GridIndex, STENCIL* Stencil, GRID_RANGE* Range)
{
	long				i = 0;
	long				j = 0;
	long				k = 0;
//	long				x = 0;
	long				y = 0;
	long				z = 0;
	short				Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Expansion >> 1;
//	double				h = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->cmn.h;
	long				iMax = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Nx - Exp - 1;
	long				jMax = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Ny - Exp - 1;
	long				kMax = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Nz - Exp - 1;
	long				X = 0;

	//	Given a grid point index, return the *valid* grid point indices corresponding to the input grid point index.
	//		-> How to handle edge cases
	//		-> Grid memory is probably more important than stencil memory
	//		-> How many slices are there?

	//	Figure out the number of "slices" (HOLD OFF - THIS SHOULD BE DONE IN THE STENCIL ITSELF)

	//	Figure out which grid point (indices?) go with each slice
	//		For z(size) -> x-y plane
	//			for y(z) -> x line
	//				Each slice is a line, return the range of x(z,y) coordinates valid for the stencil

	//	GridIndex -> (i,j,k)
	rectangular_row_major_b_spline_idx2ijk(Grid, GridIndex, &i, &j, &k);
//printf("GridIndex: %04ld (%04ld, %04ld, %04ld) --> %ld slices\n", GridIndex, i,j,k, Slices);

	Range->NumSlices = 0;
	for (z = -Stencil->Size; z <= Stencil->Size; z++)
	{
//		if (k-abs(z) < -Exp || k+abs(z) > kMax)	continue;
		if (k+z < -Exp || k+z > kMax)	continue;
		for (y = -Stencil->YMax2[abs(z)]; y <= Stencil->YMax2[abs(z)]; y++)
		{
//			if (j-abs(y) < -Exp || j+abs(y) > jMax) continue;
			if (j+y < -Exp || j+y > jMax) continue;
			X = Stencil->XMax2[STENCIL_MAP_X(MIN(abs(y),abs(z))) + STENCIL_MAP_Y(MAX(abs(z),abs(y)))];	//	Assumes y <= z
			Range->Ranges[Range->NumSlices].Min = rectangular_row_major_b_spline_ijk2idx(Grid, MAX(i-X, -Exp),j+y,k+z);
			Range->Ranges[Range->NumSlices].Max = rectangular_row_major_b_spline_ijk2idx(Grid, MIN(i+X, iMax),j+y,k+z);
			Range->NumSlices++;
		}
	}
}

void	rectangular_row_major_b_spline_get_grid_points_stencil_top(void* Grid, long GridIndex, GRID_RANGE* Range)
{
	//	Currently, the top gird calculation is "all to all"
	rectangular_row_major_b_spline_get_grid_points_all(Grid, Range);
}

double	rectangular_row_major_b_spline_get_grid_point_value(void* Grid, long GridIndex)
{
	return ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Data[GridIndex];
}

void	rectangular_row_major_b_spline_increment_grid_point_value(void* Grid, long GridIndex, double Value)
{
	((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Data[GridIndex] += Value;
}

void	rectangular_row_major_b_spline_create_copy_grid_structure(void* DstGrid, void* SrcGrid)
{
	SIMULATION_DOMAIN		Domain;
	short					Level = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) SrcGrid)->cmn.Level;
	double					h = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) SrcGrid)->cmn.h;
	long					Nx = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) SrcGrid)->Nx;
	long					Ny = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) SrcGrid)->Ny;
	long					Nz = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) SrcGrid)->Nz;
	short					Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) SrcGrid)->Expansion;

	Domain.MinimumCoordinates.x = 0.0;
	Domain.MinimumCoordinates.y = 0.0;
	Domain.MinimumCoordinates.z = 0.0;

	Domain.MaximumCoordinates.x = (double) h*(Nx - Exp - 1.0);
	Domain.MaximumCoordinates.y = (double) h*(Ny - Exp - 1.0);
	Domain.MaximumCoordinates.z = (double) h*(Nz - Exp - 1.0);

	//	Initialize just like any other grid, except Nx,Ny,Nz will be set how we want them
	grid_initialize(DstGrid, &Domain, Level, h/pow(2.0, Level), Exp);
}

void	rectangular_row_major_b_spline_create_finer_grid(void* FineGrid, void* CoarseGrid)
{
	SIMULATION_DOMAIN		Domain;
	short					Level = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) CoarseGrid)->cmn.Level;
	double					h = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) CoarseGrid)->cmn.h;
	long					Nx = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) CoarseGrid)->Nx;
	long					Ny = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) CoarseGrid)->Ny;
	long					Nz = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) CoarseGrid)->Nz;
	short					Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) CoarseGrid)->Expansion;

	Domain.MinimumCoordinates.x = 0.0;
	Domain.MinimumCoordinates.y = 0.0;
	Domain.MinimumCoordinates.z = 0.0;

	Domain.MaximumCoordinates.x = (double) h*(Nx - Exp - 1.0);
	Domain.MaximumCoordinates.y = (double) h*(Ny - Exp - 1.0);
	Domain.MaximumCoordinates.z = (double) h*(Nz - Exp - 1.0);

	//	Initialize just like any other grid, except Nx,Ny,Nz will be set how we want them
	grid_initialize(FineGrid, &Domain, Level-1, h/pow(2.0, Level), Exp);
}

void	rectangular_row_major_b_spline_create_coarser_grid(void* CoarseGrid, void* FineGrid)
{
	SIMULATION_DOMAIN		Domain;
	short					Level = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) FineGrid)->cmn.Level;
	double					h = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) FineGrid)->cmn.h;
	long					Nx = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) FineGrid)->Nx;
	long					Ny = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) FineGrid)->Ny;
	long					Nz = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) FineGrid)->Nz;
	short					Exp = ((RECTANGULAR_ROW_MAJOR_B_SPLINE*) FineGrid)->Expansion;

	Domain.MinimumCoordinates.x = 0.0;
	Domain.MinimumCoordinates.y = 0.0;
	Domain.MinimumCoordinates.z = 0.0;

	Domain.MaximumCoordinates.x = (double) h*(Nx - Exp - 1.0);
	Domain.MaximumCoordinates.y = (double) h*(Ny - Exp - 1.0);
	Domain.MaximumCoordinates.z = (double) h*(Nz - Exp - 1.0);

	//	Initialize just like any other grid, except Nx,Ny,Nz will be set how we want them
	grid_initialize(CoarseGrid, &Domain, Level+1, h/pow(2.0, Level), Exp);
}

void	rectangular_row_major_b_spline_uninitialize(void* Grid)
{
	assert(Grid != NULL);

	//	Release the grid memory
	dynfree(((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Data);
}

//	INTERNAL Methods

//	End of file
