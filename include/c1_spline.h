//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.h - C1 phi
*/

#ifndef	C1_SPLINE_H
#define	C1_SPLINE_H

#include "all.h"
#include "memory.h"
#include "interpolant.h"

typedef struct
{
	//	Members
	INTERPOLANT		cmn;
} C1_SPLINE;

//	EXTERNAL Methods
void c1_spline_initialize(void* Interpolant);
void c1_spline_evaluate(void* Interpolant);
void c1_spline_uninitialize(void* Interpolant);

//	INTERNAL Methods
void c1_spline_compute_g2p(C1_SPLINE* C1);
void c1_spline_compute_g2fg(C1_SPLINE* C1);
void c1_spline_compute_g2g(C1_SPLINE* C1);
void c1_spline_compute_tg2g(C1_SPLINE* C1);

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phiC1(int p, double x, double* dphi);

/*
Compute the coefficients for the interpolant which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of interpolant
g2fg is (p + 1)-vector
*/
void compute_g2fgC1(short p, double* g2fg);

#endif

// End of file