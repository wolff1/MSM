//-------|---------|---------|---------|---------|---------|---------|---------|
/*
grid.c - Routines for the generic GRID class
*/

#include "grid.h"

void grid_initialize(void* Grid, SIMULATION_DOMAIN* Domain, short Level, double h)
{
	assert(Grid != NULL);

	//	Initialize Members
	((GRID*) Grid)->Level = Level;
	((GRID*) Grid)->h = pow(2.0, Level)*h;

	//	Initialize Grid by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	//(*Init)(Grid, Domain);
	(*((GRID*) Grid)->initialize)(Grid, Domain);
}

void grid_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	Copy Members
	((GRID*) Dst)->h = ((GRID*) Src)->h;
	((GRID*) Dst)->Level = ((GRID*) Src)->Level;

	//	Copy Methods
	((GRID*) Dst)->initialize					= ((GRID*) Src)->initialize;
	((GRID*) Dst)->copy							= ((GRID*) Src)->copy;
	((GRID*) Dst)->xyz2idx						= ((GRID*) Src)->xyz2idx;
	((GRID*) Dst)->ijk2idx						= ((GRID*) Src)->ijk2idx;
	((GRID*) Dst)->get_grid_points_all			= ((GRID*) Src)->get_grid_points_all;
	((GRID*) Dst)->get_grid_points_coarse		= ((GRID*) Src)->get_grid_points_coarse;
	((GRID*) Dst)->get_grid_points_stencil		= ((GRID*) Src)->get_grid_points_stencil;
	((GRID*) Dst)->get_grid_points_stencil_top	= ((GRID*) Src)->get_grid_points_stencil_top;
	((GRID*) Dst)->get_grid_point_value			= ((GRID*) Src)->get_grid_point_value;
	((GRID*) Dst)->increment_grid_point_value	= ((GRID*) Src)->increment_grid_point_value;
	((GRID*) Dst)->create_finer_grid			= ((GRID*) Src)->create_finer_grid;
	((GRID*) Dst)->create_coarser_grid			= ((GRID*) Src)->create_coarser_grid;
	((GRID*) Dst)->uninitialize					= ((GRID*) Src)->uninitialize;

	//	Call sub-class copy method
	(*((GRID*) Src)->copy)(Dst, Src);
}

void grid_uninitialize(void* Grid)
{
	assert(Grid != NULL);

	((GRID*) Grid)->uninitialize(Grid);
}

//	End of file