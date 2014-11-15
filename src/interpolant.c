//-------|---------|---------|---------|---------|---------|---------|---------|
/*
interpolant.c - Parent (abstract) class for interpolants. Children must
				implement the interfaces in order for MSM to work.
*/

#include "interpolant.h"

//	EXTERNAL Methods
void interpolant_initialize(void** Interpolant, size_t Size, void* Init(void*,MSM_PARAMETERS*), MSM_PARAMETERS* MsmParams)
{	//	NOTE: Interpolant is ADDRESS of a void*
	assert(*Interpolant == NULL);
	assert(MsmParams != NULL);

	//	Dynamically allocate zero-ed out memory for *Interpolant
	*Interpolant = dynmem(Size);

	//	Initialize Members
	((INTERPOLANT*)(*Interpolant))->p = MsmParams->p;
	((INTERPOLANT*)(*Interpolant))->g2p = NULL;
	((INTERPOLANT*)(*Interpolant))->g2fg = NULL;
	stencil_initialize(&((INTERPOLANT*)(*Interpolant))->GammaI, (long) ceil(2.0*MsmParams->a / MsmParams->h), STENCIL_SHAPE_SPHERE);
	stencil_initialize(&((INTERPOLANT*)(*Interpolant))->GammaT, (long) ceil(MsmParams->D), STENCIL_SHAPE_CUBE);
	stencil_initialize(&((INTERPOLANT*)(*Interpolant))->g2g, ((INTERPOLANT*)(*Interpolant))->GammaI.size, ((INTERPOLANT*)(*Interpolant))->GammaI.shape);
	stencil_initialize(&((INTERPOLANT*)(*Interpolant))->tg2g, ((INTERPOLANT*)(*Interpolant))->GammaT.size, ((INTERPOLANT*)(*Interpolant))->GammaT.shape);

	//	Initialize Interpolant by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(*Interpolant, MsmParams);
}

//	End of file