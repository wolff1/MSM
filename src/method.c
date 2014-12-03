//-------|---------|---------|---------|---------|---------|---------|---------|
/*
method.c - Routines for the generic METHOD class
*/

#include "method.h"

void method_initialize(void* Method, void* Init(void*), short Id)
{
	assert(Method != NULL);

	//	Initialize Members
	((METHOD*) Method)->Id = Id;

	//	Initialize Method by calling function pointer to its initialize routine
	//		-> This routine MUST set the other function pointers appropriately!
	(*Init)(Method);
}

void method_copy(METHOD* SrcMethod, METHOD* DstMethod)
{
	//	Members
	DstMethod->Id = SrcMethod->Id;

	//	Methods
	DstMethod->preprocess = SrcMethod->preprocess;
	DstMethod->evaluate = SrcMethod->evaluate;
	DstMethod->uninitialize = SrcMethod->uninitialize;

	//	Call sub-class copy method
}

void method_uninitialize(void* Method)
{
	assert(Method != NULL);

	((METHOD*) Method)->uninitialize(Method);
}

//	End of file