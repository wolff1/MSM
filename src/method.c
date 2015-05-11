//-------|---------|---------|---------|---------|---------|---------|---------|
/*
method.c - Routines for the generic METHOD class
*/

#include "method.h"

void method_initialize(void* Method, void* Init(void*,void*,void*), short Id, void* Parameters, void* Options)
{
	assert(Method != NULL);

	//	Initialize Members
	((METHOD*) Method)->Id = Id;
	((METHOD*) Method)->Size = sizeof(METHOD);

	//	Initialize Method by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(Method, Parameters, Options);
}

void method_copy(METHOD* Dst, METHOD* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	Copy Members
	Dst->Id = Src->Id;
	Dst->Size = Src->Size;

	//	Copy Methods
	Dst->copy = Src->copy;
	Dst->preprocess = Src->preprocess;
	Dst->evaluate = Src->evaluate;
	Dst->uninitialize = Src->uninitialize;

	//	Call sub-class copy method
	(*Src->copy)((void*)Dst, (void*)Src);
}

void method_uninitialize(void* Method)
{
	assert(Method != NULL);

	((METHOD*) Method)->uninitialize(Method);
}

//	End of file