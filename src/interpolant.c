//-------|---------|---------|---------|---------|---------|---------|---------|
/*
interpolant.c - Parent (abstract) class for interpolants. Children must
				implement the interfaces in order for MSM to work.
*/

#include "interpolant.h"

//	EXTERNAL Methods
void interpolant_initialize(void* Interpolant, void* Init(void*,MSM_PARAMETERS*), MSM_PARAMETERS* MsmParams)
{	//	NOTE: Interpolant is ADDRESS of a void*
	assert(Interpolant != NULL);
	assert(MsmParams != NULL);

	//	Initialize Members
	((INTERPOLANT*)Interpolant)->p = MsmParams->p;
	((INTERPOLANT*)Interpolant)->g2p = NULL;
	((INTERPOLANT*)Interpolant)->g2fg = NULL;
	((INTERPOLANT*)Interpolant)->g2g = NULL;
	((INTERPOLANT*)Interpolant)->tg2g = NULL;
	((INTERPOLANT*)Interpolant)->Size = sizeof(INTERPOLANT);

	//	Initialize Interpolant by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(Interpolant, MsmParams);
}

void interpolant_copy(INTERPOLANT* Dst, INTERPOLANT* Src)
{
	short		i = 0;

	assert(Dst != NULL);
	assert(Src != NULL);

	//	Copy Members
	Dst->p = Src->p;

	Dst->g2p = (double**) dynarr_d(Dst->p/2, Dst->p);
	for (i = 0; i < Dst->p/2; i++)
	{
		memcpy(Dst->g2p[i], Src->g2p[i], (Dst->p)*sizeof(double));
	}

	Dst->g2fg = (double*) dynvec(Dst->p/2+1, sizeof(double));
	memcpy(Dst->g2fg, Src->g2fg, (Dst->p/2+1)*sizeof(double));

	Dst->g2g = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_copy(Dst->g2g, Src->g2g);

	Dst->tg2g = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_copy(Dst->tg2g, Src->tg2g);

	Dst->Size = Src->Size;

	//	Copy Methods
	Dst->copy = Src->copy;
	Dst->compute_g2g = Src->compute_g2g;
	Dst->compute_tg2g = Src->compute_tg2g;
	Dst->evaluate = Src->evaluate;
	Dst->uninitialize = Src->uninitialize;

	//	Call sub-class copy method
	(*Src->copy)((void*)Dst, (void*)Src);
}

//	End of file