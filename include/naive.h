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
#include "particle.h"

typedef struct
{
	//	Members
	METHOD				cmn;
} NAIVE;

//	EXTERNAL Methods
void naive_initialize(void* Naive);
void naive_copy(void* Dst, void* Src);
void naive_preprocess(void* Naive, double DomainRadius);
void naive_evaluate(void* Naive, long N, PARTICLE* r, double* U, double** f);
void naive_uninitialize(void* Naive);

//	INTERNAL Methods

#endif

//	End of file