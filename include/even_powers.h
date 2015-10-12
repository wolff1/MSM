//-------|---------|---------|---------|---------|---------|---------|---------|
/*
even_powers.h - the softening function used to split the kernel
*/

#ifndef	EVEN_POWERS_H
#define	EVEN_POWERS_H

#include "all.h"
#include "softener.h"
#include "memory.h"

typedef struct
{
	SOFTENER		cmn;
} EVEN_POWERS;

//	EXTERNAL Methods
void even_powers_initialize(void* Softener);
void even_powers_copy(void* Dst, void* Src);
void even_powers_soften(void* Softener, long Len, double* X, double* F, double* DF);
void even_powers_split(void* Softener, long Len, double* X, double* F, double* DF);
void even_powers_derivative(void* Softener, long Len, double* X, double* F, short DerivativeNumber);
void even_powers_uninitialize(void* Softener);

//	INTERNAL Methods
void even_powers_compute_p2p(EVEN_POWERS* Ep);

#endif

// End of file