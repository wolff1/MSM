//-------|---------|---------|---------|---------|---------|---------|---------|
/*
softener.c - 
*/

#include "softener.h"

void softener_initialize(void** Softener, size_t Size, void* Init(void*), short k)
{	//	NOTE: Softener is ADDRESS of a void*
	assert(*Softener == NULL);

	//	Dynamically allocate zero-ed out memory for *Softener
	*Softener = dynmem(Size);

	//	Initialize Members
	((SOFTENER*)(*Softener))->k = k;
	((SOFTENER*)(*Softener))->p2p = NULL;

	//	Initialize Softner by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(*Softener);
}

//	End of file