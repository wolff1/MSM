//-------|---------|---------|---------|---------|---------|---------|---------|
/*
method.c - Routines for the generic METHOD class
*/

#include "method.h"

void method_initialize(void** Method, size_t Size, void* Init(void*))
{	//	NOTE: Method is ADDRESS of a void*
	assert(*Method == NULL);

	//	Dynamically allocate zero-ed out memory for *Method
	*Method = dynmem(Size);

	//	Initialize Members
	((METHOD*) *Method)->U = 0.0;
	((METHOD*) *Method)->f = NULL;	//	Allocated in [method]_evaluate()

	//	Initialize Method by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(*Method);
}

//	End of file