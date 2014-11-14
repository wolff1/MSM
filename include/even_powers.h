//-------|---------|---------|---------|---------|---------|---------|---------|
/*
even_powers.h - the softening function used to split the kernel
*/

#include "all.h"
#include "memory.h"
#include "output.h"			//	remove?
#include "interpolant.h"	//	remove?
#include "b_spline.h"		//	remove?

#define GAMMA_DATA			"gamma.dat"
#define THETA_DATA			"theta.dat"

#define GAMMA_DATA_LEN		strlen(GAMMA_DATA)
#define THETA_DATA_LEN		strlen(THETA_DATA)

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

// End of file