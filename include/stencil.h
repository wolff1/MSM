//-------|---------|---------|---------|---------|---------|---------|---------|
/*
stencil.h - 
*/

#ifndef STENCIL_H
#define	STENCIL_H

#include "all.h"
#include "memory.h"
#include "softener.h"

#define	STENCIL_STORAGE(L)		(L+1)*(L+2)*(L+3)/6
#define	STENCIL_MAP_Z(z)		(z)*(z+1)*(z+2)/6
#define	STENCIL_MAP_Y(y)		(y)*(y+1)/2
#define STENCIL_MAP_X(x)		(x)

//	The following are for the "stacked" symmetric plane stencils
//	which are intermediate results while applying the anti-blurring
//	operator.
#define	STENCIL_MAP_Z2(L,z)		(z)*(L+2)*(L+1)/2	//	L is K->size
#define	STENCIL_MAP_Y2(y)		(y)*(y+1)/2
#define STENCIL_MAP_X2(x)		(x)

#define	STENCIL_STORAGE_2D(L)	(L+1)*(L+2)/2	//	Used for loop ranges

#define	STENCIL_SHAPE_SPHERE	1
#define	STENCIL_SHAPE_CUBE		2

#define	STENCIL_FUNCTION_TYPE_THETA		1
#define	STENCIL_FUNCTION_TYPE_GAMMA		2

typedef struct
{
	short			Shape;
	long			Size;	//	= zmax
	long*			YMax;
	long*			XMax;
	double*			Data;
} STENCIL;

void stencil_initialize(STENCIL** t, long Size, short Shape);
void stencil_populate(STENCIL* s, SOFTENER* Softener, short FunctionType, double Scale);
void stencil_shift(STENCIL* s, short Degree, double* OmegaPrime, STENCIL* K);
void stencil_shift_infinite(STENCIL* s, short Degree, double* OmegaPrime, STENCIL* K);
void stencil_display(STENCIL* s, double h_a);
void stencil_free(STENCIL* s);

#endif

//	End of file