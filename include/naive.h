//-------|---------|---------|---------|---------|---------|---------|---------|
/*
naive.h -
*/

#ifndef	NAIVE_H
#define	NAIVE_H

#include "all.h"
#include "interpolant.h"
#include "softener.h"
#include "method.h"

typedef struct
{
	//	Members
	METHOD				cmn;
} NAIVE;

//	EXTERNAL Methods
void naive_initialize(void* Naive);
void naive_preprocess(void* Naive);
void naive_evaluate(void* Naive);
void naive_uninitialize(void* Naive);

//	INTERNAL Methods

#endif

//	End of file