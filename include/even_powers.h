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
void even_powers_uninitialize(void* Softener);

//	INTERNAL Methods
/*
Solves for coefficients of gamma
k is degree of continuity
x is k+1 vector to hold coefficients
*/
void gamma_init(short k, double* x);

/*
Evaluate gamma and gamma' at position x
c is coefficient vector for gamma
*/
double _gamma(double *c, short k, double x, double* dgamma);

/*
Evaluate theta* and theta*' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta_star(double *c, short k, double x, double* dtheta_star);

/*
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta(double *c, short k, double x, double* dtheta);

/*
Evaluate theta and theta' at position x as single polynomial
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double thetap(double *c, short k, double x, double* dtheta);

#endif

// End of file