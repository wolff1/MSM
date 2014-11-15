//-------|---------|---------|---------|---------|---------|---------|---------|
/*
softener.c - 
*/

#include "softener.h"

void softener_initialize(void* Softener, size_t Size, void* Init(void*))
{	
	assert(Softener != NULL);

	//	Initialize Softner by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(Softener);
}

//	End of file