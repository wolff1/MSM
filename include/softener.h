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
	void		(*gamma);
	void		(*theta);
} SOFTENER;

//	EXTERNAL Methods
void softener_initialize(void* Softener, size_t Size, void* Init(void*));

#endif

//	End of file