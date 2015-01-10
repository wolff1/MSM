//-------|---------|---------|---------|---------|---------|---------|---------|
/*
grid.h - 
*/

#ifndef	GRID_H
#define	GRID_H

#include "all.h"

typedef struct
{
	long		Min;
	long		Max;
} GRID_RANGE_MIN_MAX;

typedef struct
{
	long					NumSlices;
	GRID_RANGE_MIN_MAX*		Ranges;
} GRID_RANGE;

typedef struct
{
	//	Members
	double*		data;
	double		h;
	short		level;
	long		Nx;
	long		Ny;
	long		Nz;
	short		NumIntergridCoefficients;

	//	Methods
	long	(*xyz2idx)(void*, double, double, double);
	long	(*ijk2idx)(void*, long, long, long);
	void	(*get_grid_points_all)(void*, GRID_RANGE*);
	void	(*get_grid_points_coarse)(void*, long, GRID_RANGE*);
	void	(*get_grid_points_stencil)(void*, long, GRID_RANGE*);
	void	(*get_grid_points_stencil_top)(void*, long, GRID_RANGE*);
	double	(*get_grid_point_value)(void*, long);
	void	(*increment_grid_point_value)(void*, long, double value);
} GRID;

//	EXTERNAL Methods
void grid_initialize(void* Grid, void* Init(void));
void grid_copy(GRID* Dst, GRID* Src);

//	INTERNAL Methods

//	INTERNAL HELPER Methods

//	INTERNAL TESTING Methods

#endif