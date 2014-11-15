//-------|---------|---------|---------|---------|---------|---------|---------|
/*
centered b-spline phi
*/

#include "interpolant.h"

//	EXTERNAL Methods
void interpolant_initialize(void** Interpolant, size_t Size, void* Init(void*), short p)
{	//	NOTE: Interpolant is ADDRESS of a void*
	assert(*Interpolant == NULL);

	//	Dynamically allocate zero-ed out memory for *Interpolant
	assert(*Interpolant == NULL);

	*Interpolant = dynmem(Size);

	//	Initialize Members
	((INTERPOLANT*)(*Interpolant))->p = p;
	((INTERPOLANT*)(*Interpolant))->g2p = NULL;
	((INTERPOLANT*)(*Interpolant))->g2fg = NULL;
	((INTERPOLANT*)(*Interpolant))->g2g = NULL;
	((INTERPOLANT*)(*Interpolant))->tg2g = NULL;

	//	Initialize Interpolant by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(*Interpolant);
}

//	End of file