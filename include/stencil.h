//-------|---------|---------|---------|---------|---------|---------|---------|
/*
stencil.h - 
*/

#include "all.h"
#include "memory.h"
#include "even_powers.h"

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
	short			shape;
	long			size;	//	= zmax
	long*			ymax;
	long*			xmax;
	double*			data;
} STENCIL;

/*
stencil operations below:
*/
void stencil_initialize(STENCIL* s, long size, short shape);
void stencil_populate(STENCIL* s, double* c, short k, short function_type,
						double h_a);
void stencil_display(STENCIL* s, double h_a);
void stencil_shift(STENCIL* s, short degree, double* omega_prime, STENCIL* K);
void stencil_free(STENCIL* s);
void stencil_naive(short p, double a, double h, short degree, double* omegap, short k, double* c, STENCIL* Ki);
void stencil_shift_infinite(STENCIL* s, short degree, double* omegap, STENCIL* K);
void stencil_naive_top(short p, double a, double h, short degree, double* omegap, short k, double* c, STENCIL* Kt);

//	End of file