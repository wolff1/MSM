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

//	printf("Min(%4.2f,%4.2f,%4.2f) Max(%4.2f,%4.2f,%4.2f)\n",
//		Domain->MinimumCoordinates.x,Domain->MinimumCoordinates.y,Domain->MinimumCoordinates.z,
//		Domain->MaximumCoordinates.x,Domain->MaximumCoordinates.y,Domain->MaximumCoordinates.z);

	//	Initialize COMMON members

	//	Initialize COMMON function pointers
	MyGrid->cmn.copy							= &rectangular_row_major_b_spline_copy;
	MyGrid->cmn.display							= &rectangular_row_major_b_spline_display;
	MyGrid->cmn.xyz2idx							= &rectangular_row_major_b_spline_xyz2idx;
	MyGrid->cmn.ijk2idx							= &rectangular_row_major_b_spline_ijk2idx;
	MyGrid->cmn.get_grid_points_all_max_slices	= &rectangular_row_major_b_spline_get_grid_points_all_max_slices;
	MyGrid->cmn.get_grid_points_all				= &rectangular_row_major_b_spline_get_grid_points_all;
	MyGrid->cmn.get_grid_points_coarse			= &rectangular_row_major_b_spline_get_grid_points_coarse;
	MyGrid->cmn.get_grid_points_stencil			= &rectangular_row_major_b_spline_get_grid_points_stencil;
	MyGrid->cmn.get_grid_points_stencil_top		= &rectangular_row_major_b_spline_get_grid_points_stencil_top;
	MyGrid->cmn.get_grid_point_value			= &rectangular_row_major_b_spline_get_grid_point_value;
	MyGrid->cmn.increment_grid_point_value		= &rectangular_row_major_b_spline_increment_grid_point_value;
	MyGrid->cmn.create_finer_grid				= &rectangular_row_major_b_spline_create_finer_grid;
	MyGrid->cmn.create_coarser_grid				= &rectangular_row_major_b_spline_create_coarser_grid;
	MyGrid->cmn.uninitialize					= &rectangular_row_major_b_spline_uninitialize;

	//	Initialize RECTANGULAR_ROW_MAJOR_B_SPLINE members
	MyGrid->Expansion = Expansion;
	MyGrid->Nx = (long) 2*ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / (2*MyGrid->cmn.h)) + Expansion + 1;
	MyGrid->Ny = (long) 2*ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / (2*MyGrid->cmn.h)) + Expansion + 1;
	MyGrid->Nz = (long) 2*ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / (2*MyGrid->cmn.h)) + Expansion + 1;

	MyGrid->Data = (double*) dynvec(MyGrid->Nx*MyGrid->Ny*MyGrid->Nz, sizeof(double));

	MyGrid->NumIntergridCoefficients = Expansion/2 + 1;
}

void	rectangular_row_major_b_spline_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	COMMON items copied in grid_copy
//	double*		Data;
//	long		Nx;
//	long		Ny;
//	long		Nz;
//	short		NumIntergridCoefficients;
}

void	rectangular_row_major_b_spline_display(void* Grid)
{
	long								i = 0;
	long								j = 0;
	long								k = 0;
	RECTANGULAR_ROW_MAJOR_B_SPLINE*		MyGrid = (RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid;

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

void	rectangular_row_major_b_spline_get_grid_points_coarse(void* Grid, long FineGridIndex, GRID_RANGE* Range)
{
}

void	rectangular_row_major_b_spline_get_grid_points_stencil(void* Grid, long GridIndex, GRID_RANGE* Range)
{
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

void	rectangular_row_major_b_spline_create_finer_grid(void* FineGrid, void* CoarseGrid)
{
}

void	rectangular_row_major_b_spline_create_coarser_grid(void* CoarseGrid, void* FineGrid)
{
}

void	rectangular_row_major_b_spline_uninitialize(void* Grid)
{
	assert(Grid != NULL);

	//	Release the grid memory
	dynfree(((RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid)->Data);
}

//	INTERNAL Methods

//	End of file
