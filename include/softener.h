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
	void		(*soften)(void*, long, double*, double*, double*);
	void		(*split)(void*, long, double*, double*, double*);
	void		(*uninitialize)(void*);
} SOFTENER;

//	EXTERNAL Methods
void softener_initialize(void* Softener, void* Init(void*), short k);

#endif

//	End of file