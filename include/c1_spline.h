//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.h - C1 phi
*/

#include "all.h"
#include "memory.h"
#include "output.h"		//	remove?

#define PHI_DATA_C1		"phiC1.dat"

#define PHI_DATA_C1_LEN	strlen(PHI_DATA_C1)

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

// End of file