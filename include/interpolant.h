//-------|---------|---------|---------|---------|---------|---------|---------|
/*
interpolant.h - centered b-spline phi
*/

#include "all.h"
#include "memory.h"
#include "output.h"			//	remove?

#define PHI_DATA		"phi.dat"

#define PHI_DATA_LEN	strlen(PHI_DATA)

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(short p, double x, double* dphi);

void new_phi(short p, double** g2p, short n, double* x, double* phi, double* dphi);

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of B-spline
g2fg is (p/2 + 1)-vector
*/
void compute_g2fg(short p, double* g2fg);

// End of file