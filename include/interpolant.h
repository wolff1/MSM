//-------|---------|---------|---------|---------|---------|---------|---------|
/*
interpolant.h - centered b-spline phi
*/

#ifndef	INTERPOLANT_H
#define	INTERPOLANT_H

#include "all.h"
#include "memory.h"

typedef struct
{
	//	Members
	short		p;
	double*		g2p;
	double*		g2fg;
	double*		g2g;
	double*		tg2g;

	//	Methods
	void		(*evaluate)(void*);
	void		(*uninitialize)(void*);
} INTERPOLANT;

//	EXTERNAL Methods
void interpolant_initialize(void** Interpolant, size_t Size, void* Init(void*), short p);

#endif

// End of file