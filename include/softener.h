//-------|---------|---------|---------|---------|---------|---------|---------|
/*
softener.h - 
*/

#ifndef	SOFTENER_H
#define	SOFTENER_H

#include "all.h"

typedef struct
{
	//	Members
	short		k;
	double*		p2p;
	//	Methods
	void*		initialize;
	void*		gamma;
	void*		theta;
} SOFTENER;

#endif

//	End of file