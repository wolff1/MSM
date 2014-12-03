//-------|---------|---------|---------|---------|---------|---------|---------|
/*
method.h - Definition for the generic METHOD class
*/

#ifndef	METHOD_H
#define	METHOD_H

#include "all.h"
#include "memory.h"

typedef struct
{
	//	Members
	short		Id;

	//	Methods
	void		(*preprocess)(void*,double);
	void		(*evaluate)(void*);
	void		(*uninitialize)(void*);
} METHOD;

//	EXTERNAL Methods
void method_initialize(void* Method, void* Init(void*), short Id);
void method_copy(METHOD* SrcMethod, METHOD* DstMethod);
void method_uninitialize(void* Method);

#endif

//	End of file