//-------|---------|---------|---------|---------|---------|---------|---------|
/*
rectangular_row_major_b_spline.h - 
*/

#ifndef	RECTANGULAR_ROW_MAJOR_B_SPLINE_H
#define	RECTANGULAR_ROW_MAJOR_B_SPLINE_H

#include "all.h"
#include "grid.h"

typedef struct
{
	//	COMMON Members
	GRID		cmn;

	//	MY Members
	long		Nx;
	long		Ny;
	long		Nz;
	short		NumIntergridCoefficients;
	short		Expansion;
	double*		Data;
} RECTANGULAR_ROW_MAJOR_B_SPLINE;

//	EXTERNAL Methods
void	rectangular_row_major_b_spline_initialize(void* Grid, SIMULATION_DOMAIN* Domain, short Expansion);
void	rectangular_row_major_b_spline_copy(void* DstGrid, void* SrcGrid);
void	rectangular_row_major_b_spline_display(void* Grid);
long	rectangular_row_major_b_spline_xyz2idx(void* Grid, SIMULATION_DOMAIN* Domain, double x, double y , double z);
long	rectangular_row_major_b_spline_ijk2idx(void* Grid, long i, long j, long k);
void	rectangular_row_major_b_spline_idx2ijk(void* Grid, long Idx, long* i, long* j, long* k);
long	rectangular_row_major_b_spline_get_grid_points_all_max_slices(void* Grid);
void	rectangular_row_major_b_spline_get_grid_points_all(void* Grid, GRID_RANGE* Range);
void	rectangular_row_major_b_spline_get_grid_points_coarse(void* Grid, long FineGridIndex, GRID_RANGE* Range);
void	rectangular_row_major_b_spline_get_grid_points_stencil(void* Grid, long GridIndex, STENCIL* Stencil, GRID_RANGE* Range);
void	rectangular_row_major_b_spline_get_grid_points_stencil_top(void* Grid, long GridIndex, GRID_RANGE* Range);
double	rectangular_row_major_b_spline_get_grid_point_value(void* Grid, long GridIndex);
void	rectangular_row_major_b_spline_increment_grid_point_value(void* Grid, long GridIndex, double Value);
void	rectangular_row_major_b_spline_create_copy_grid_structure(void* DstGrid, void* SrcGrid);
void	rectangular_row_major_b_spline_create_finer_grid(void* FineGrid, void* CoarseGrid);
void	rectangular_row_major_b_spline_create_coarser_grid(void* CoarseGrid, void* FineGrid);
void	rectangular_row_major_b_spline_uninitialize(void* Grid);

//	INTERNAL Methods

//	INTERNAL HELPER Methods

//	INTERNAL TESTING Methods

#endif