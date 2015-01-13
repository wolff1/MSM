//-------|---------|---------|---------|---------|---------|---------|---------|
/*
grid.h - 
*/

#ifndef	GRID_H
#define	GRID_H

#include "all.h"
#include "memory.h"
#include "simulation_domain.h"

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
	short		Level;
	double		h;		//	Actually h_l

	//	Methods
	void		(*initialize)(void*, SIMULATION_DOMAIN*, short);
	void		(*copy)(void*, void*);
	long		(*xyz2idx)(void*, SIMULATION_DOMAIN*, double, double, double);
	long		(*ijk2idx)(void*, long, long, long);
	void		(*get_grid_points_all)(void*, GRID_RANGE*);
	void		(*get_grid_points_coarse)(void*, long, GRID_RANGE*);
	void		(*get_grid_points_stencil)(void*, long, GRID_RANGE*);
	void		(*get_grid_points_stencil_top)(void*, long, GRID_RANGE*);
	double		(*get_grid_point_value)(void*, long);
	void		(*increment_grid_point_value)(void*, long, double);
	void		(*create_finer_grid)(void*, void*);
	void		(*create_coarser_grid)(void*, void*);
	void		(*uninitialize)(void*);
} GRID;

//	EXTERNAL Methods
void grid_initialize(void* Grid, SIMULATION_DOMAIN* Domain, short Level, double h, short Expansion);
void grid_copy(void* Dst, void* Src);
void grid_uninitialize(void* Grid);

//	INTERNAL Methods

//	INTERNAL HELPER Methods

//	INTERNAL TESTING Methods

#endif