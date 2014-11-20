//-------|---------|---------|---------|---------|---------|---------|---------|
/*
softener.c - 
*/

#include "softener.h"

void softener_initialize(void* Softener, void* Init(void*), short k)
{
	assert(Softener != NULL);

	//	Initialize Members
	((SOFTENER*) Softener)->k = k;
	((SOFTENER*) Softener)->p2p = NULL;

	//	Initialize Softner by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(Softener);
}

//	End of file