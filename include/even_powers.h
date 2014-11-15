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
void even_powers_soften(void* Softener, long Len, double* X, double* F, double* DF);
void even_powers_split(void* Softener, long Len, double* X, double* F, double* DF);
void even_powers_uninitialize(void* Softener);

//	INTERNAL Methods
void even_powers_compute_p2p(EVEN_POWERS* Ep);

//*****************************************************//
//*****DELETE EVERYTHING FROM HERE DOWN EVENTUALLY*****//
//*****************************************************//

void gamma_init(short k, double* x);
double _gamma(double *c, short k, double x, double* dgamma);
double theta_star(double *c, short k, double x, double* dtheta_star);
double theta(double *c, short k, double x, double* dtheta);
double thetap(double *c, short k, double x, double* dtheta);

#endif

// End of file