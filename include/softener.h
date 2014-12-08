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
	size_t		Size;

	//	Methods
	void		(*copy)(void*,void*);
	void		(*soften)(void*, long, double*, double*, double*);
	void		(*split)(void*, long, double*, double*, double*);
	void		(*uninitialize)(void*);
} SOFTENER;

//	EXTERNAL Methods
void softener_initialize(void* Softener, void* Init(void*), short k);
void softener_copy(void* Dst, void* Src);

#endif

//	End of file