//-------|---------|---------|---------|---------|---------|---------|---------|
/*
softener.h - 
*/

#ifndef	SOFTENER_H
#define	SOFTENER_H

#include "all.h"
#include "memory.h"

typedef struct
{
	//	Members
	short		k;
	double*		p2p;

	//	Methods
	double		(*gamma)(double*, short, double, double*);
	double		(*theta)(double*, short, double, double*);
	void		(*uninitialize)(void*);
} SOFTENER;

//	EXTERNAL Methods
void softener_initialize(void** Softener, size_t Size, void* Init(void*), short k);

#endif

//	End of file