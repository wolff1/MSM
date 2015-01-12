//-------|---------|---------|---------|---------|---------|---------|---------|
/*
rectangular_row_major_b_spline.c -
*/

#include "rectangular_row_major_b_spline.h"

//	EXTERNAL Methods
void	rectangular_row_major_b_spline_initialize(void* Grid, SIMULATION_DOMAIN* Domain /*, Phi_Len, alpha, stencil_shape*/)
{
	RECTANGULAR_ROW_MAJOR_B_SPLINE*		MyGrid = (RECTANGULAR_ROW_MAJOR_B_SPLINE*) Grid;
	short								p = 4;//Msm->prm.p;

	assert(MyGrid != NULL);

	//	Initialize COMMON members

	//	Initialize COMMON function pointers
	MyGrid->cmn.copy						= &rectangular_row_major_b_spline_copy;
	MyGrid->cmn.xyz2idx						= &rectangular_row_major_b_spline_xyz2idx;
	MyGrid->cmn.ijk2idx						= &rectangular_row_major_b_spline_ijk2idx;
	MyGrid->cmn.get_grid_points_all			= &rectangular_row_major_b_spline_get_grid_points_all;
	MyGrid->cmn.get_grid_points_coarse		= &rectangular_row_major_b_spline_get_grid_points_coarse;
	MyGrid->cmn.get_grid_points_stencil		= &rectangular_row_major_b_spline_get_grid_points_stencil;
	MyGrid->cmn.get_grid_points_stencil_top	= &rectangular_row_major_b_spline_get_grid_points_stencil_top;
	MyGrid->cmn.get_grid_point_value		= &rectangular_row_major_b_spline_get_grid_point_value;
	MyGrid->cmn.increment_grid_point_value	= &rectangular_row_major_b_spline_increment_grid_point_value;
	MyGrid->cmn.create_finer_grid			= &rectangular_row_major_b_spline_create_finer_grid;
	MyGrid->cmn.create_coarser_grid			= &rectangular_row_major_b_spline_create_coarser_grid;
	MyGrid->cmn.uninitialize				= &rectangular_row_major_b_spline_uninitialize;

	//	Initialize RECTANGULAR_ROW_MAJOR_B_SPLINE members
	//		-> Need:
	//			[X]	simulation domain:	for grid shape/size
	//			[ ]	stencil info:		for stencil shape/size (cube/sphere, 2*alpha)
	//			[X]	interpolant:		for "phi" stencil shape/size (p)

	MyGrid->Nx = (long) 2*ceil((Domain->MaximumCoordinates.x - Domain->MinimumCoordinates.x) / (2*MyGrid->cmn.h)) + p + 1;
	MyGrid->Ny = (long) 2*ceil((Domain->MaximumCoordinates.y - Domain->MinimumCoordinates.y) / (2*MyGrid->cmn.h)) + p + 1;
	MyGrid->Nz = (long) 2*ceil((Domain->MaximumCoordinates.z - Domain->MinimumCoordinates.z) / (2*MyGrid->cmn.h)) + p + 1;
printf("Grid %hd: (%ld,%ld,%ld) -> %ld\n", MyGrid->cmn.Level, MyGrid->Nx, MyGrid->Ny, MyGrid->Nz, MyGrid->Nx*MyGrid->Ny*MyGrid->Nz);
	MyGrid->Data = (double*) dynvec(MyGrid->Nx*MyGrid->Ny*MyGrid->Nz, sizeof(double));
	MyGrid->NumIntergridCoefficients = p;
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

long	rectangular_row_major_b_spline_xyz2idx(void* Grid, double x, double y , double z)
{
	return 0;
}

long	rectangular_row_major_b_spline_ijk2idx(void* Grid, long i, long j, long k)
{
	return 0;
}

void	rectangular_row_major_b_spline_get_grid_points_all(void* Grid, GRID_RANGE* Range)
{
}

void	rectangular_row_major_b_spline_get_grid_points_coarse(void* Grid, long FineGridIndex, GRID_RANGE* Range)
{
}

void	rectangular_row_major_b_spline_get_grid_points_stencil(void* Grid, long GridIndex, GRID_RANGE* Range)
{
}

void	rectangular_row_major_b_spline_get_grid_points_stencil_top(void* Grid, long GridIndex, GRID_RANGE* Range)
{
}

double	rectangular_row_major_b_spline_get_grid_point_value(void* Grid, long GridIndex)
{
	return 0.0;
}

void	rectangular_row_major_b_spline_increment_grid_point_value(void* Grid, long GridIndex, double Value)
{
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
