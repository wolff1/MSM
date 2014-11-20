//-------|---------|---------|---------|---------|---------|---------|---------|
/*
method.c - Routines for the generic METHOD class
*/

#include "method.h"

void method_initialize(void** Method, size_t Size, void* Init(void*), short Id)
{	//	NOTE: Method is ADDRESS of a void*
	assert(*Method == NULL);

	//	Dynamically allocate zero-ed out memory for *Method
	*Method = dynmem(Size);

	//	Initialize Members
	((METHOD*) *Method)->Id = Id;
	((METHOD*) *Method)->U = 0.0;
	((METHOD*) *Method)->f = NULL;	//	Allocated in [method]_evaluate()

	//	Initialize Method by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(*Method);
}

void method_copy(METHOD* SrcMethod, METHOD* DstMethod)
{
	//	Members
	DstMethod->Id = SrcMethod->Id;
	DstMethod->U = SrcMethod->U;
	DstMethod->f = SrcMethod->f;
	//	Methods
	DstMethod->preprocess = SrcMethod->preprocess;
	DstMethod->evaluate = SrcMethod->evaluate;
	DstMethod->uninitialize = SrcMethod->uninitialize;
}

void method_uninitialize(void* Method)
{
	assert(Method != NULL);

	((METHOD*) Method)->uninitialize(Method);	//	<--- FIXME - THIS MAY NOT BE A GOOD IDEA
}

//	End of file