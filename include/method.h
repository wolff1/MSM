//-------|---------|---------|---------|---------|---------|---------|---------|
/*
method.h - Definition for the generic METHOD class
*/

#ifndef	METHOD_H
#define	METHOD_H

#include "all.h"
#include "memory.h"
#include "simulation_domain.h"
//#include "particle.h"

typedef struct
{
	//	Members
	short		Id;
	size_t		Size;

	//	Methods
	void		(*copy)(void*,void*);
	void		(*preprocess)(void*,double);
	void		(*evaluate)(void*, SIMULATION_DOMAIN*);
	void		(*uninitialize)(void*);
} METHOD;

//	EXTERNAL Methods
void method_initialize(void* Method, void* Init(void*), short Id);
void method_copy(METHOD* Dst, METHOD* Src);
void method_uninitialize(void* Method);

#endif

//	End of file