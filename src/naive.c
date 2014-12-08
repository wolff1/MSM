//-------|---------|---------|---------|---------|---------|---------|---------|
/*
naive.c -
*/

#include "naive.h"

//	EXTERNAL Methods
void naive_initialize(void* Method)
{
	NAIVE*		Naive = (NAIVE*) Method;

	assert(Naive != NULL);
	printf("Initializing NAIVE!\n");

	//	Initialize COMMON members
	Naive->cmn.Size = sizeof(NAIVE);

	//	Initialize COMMON function pointers
	Naive->cmn.copy = &naive_copy;
	Naive->cmn.preprocess = &naive_preprocess;
	Naive->cmn.evaluate = &naive_evaluate;
	Naive->cmn.uninitialize = &naive_uninitialize;
}

void naive_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	--> METHOD is copied in method_copy()
}

void naive_preprocess(void* Method, double DomainRadius)
{
	NAIVE*		Naive = (NAIVE*) Method;
	assert(Naive != NULL);
	printf("NAIVE Preprocessing!\n");
}

void naive_evaluate(void* Method)
{
	NAIVE*		Naive = (NAIVE*) Method;
	assert(Naive != NULL);
	printf("NAIVE Evaluation!\n");
}

void naive_uninitialize(void* Method)
{
	NAIVE*		Naive = (NAIVE*) Method;
	assert(Naive != NULL);
	printf("Un-initializing NAIVE!\n");
}

//	INTERNAL Methods

//	End of file
