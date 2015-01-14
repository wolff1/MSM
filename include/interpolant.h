//-------|---------|---------|---------|---------|---------|---------|---------|
/*
interpolant.h - Parent (abstract) class for interpolants. Children must
				implement the interfaces in order for MSM to work.
*/

#ifndef	INTERPOLANT_H
#define	INTERPOLANT_H

#include "all.h"
#include "memory.h"
#include "stencil.h"

//	THERE MUST BE A BETTER PLACE FOR THIS THAN HERE
//		-> Conflict because needed in MSM and INTERPOLANT
typedef struct
{
	double		a;			//	Cut-off distance
	double		h;			//	Finest grid spacing
	double		alpha;		//	a/h
	short		p;			//	Gives order of interpolation accuracy
	short		k;			//	Gives degree of continuity for softening function
	short		mu;			//	Gives number of terms in series expansion for B-spline interpolant
	double		D;			//	Size of the domain
	short		L;			//	Number of *grids* used in the method
} MSM_PARAMETERS;

typedef struct
{
	//	Members
	short		p;
	double**	g2p;
	double*		g2fg;
	STENCIL*	g2g;
	STENCIL*	tg2g;
	size_t		Size;

	//	Methods
	void		(*copy)(void*,void*);
	void		(*compute_g2g)(void*, SOFTENER*, MSM_PARAMETERS*);
	void		(*compute_tg2g)(void*, SOFTENER*, MSM_PARAMETERS*);
	void		(*evaluate)(void*, long, double*, double*, double*);
	void		(*uninitialize)(void*);
} INTERPOLANT;

//	EXTERNAL Methods
void interpolant_initialize(void* Interpolant, void* Init(void*,MSM_PARAMETERS*), MSM_PARAMETERS* MsmParams);
void interpolant_copy(INTERPOLANT* Dst, INTERPOLANT* Src);

#endif

// End of file