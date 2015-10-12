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
	((SOFTENER*) Softener)->Size = sizeof(SOFTENER);

	//	Initialize Softner by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(Softener);
}

void softener_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	Copy Members
	((SOFTENER*)Dst)->k = ((SOFTENER*)Src)->k;
	//	NOTE: copy p2p in sub-class b/c its length may be dependent on the type of softener
	((SOFTENER*)Dst)->Size = ((SOFTENER*)Src)->Size;

	//	Copy Methods
	((SOFTENER*)Dst)->copy = ((SOFTENER*)Src)->copy;
	((SOFTENER*)Dst)->soften = ((SOFTENER*)Src)->soften;
	((SOFTENER*)Dst)->split = ((SOFTENER*)Src)->split;
	((SOFTENER*)Dst)->derivative = ((SOFTENER*)Src)->derivative;
	((SOFTENER*)Dst)->uninitialize = ((SOFTENER*)Src)->uninitialize;

	//	Call sub-class copy method
	(*((SOFTENER*)Src)->copy)((void*)Dst, (void*)Src);
}

//	End of file