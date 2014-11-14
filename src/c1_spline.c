//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.c - C1 phi (CUBIC ONLY)
*/

#include "c1_spline.h"

// These are for the CUBIC "shifted-powers" version
double c[2][4] = {{0.0, -0.5, 2.0, 1.5},
				  {0.0, 0.0, -0.5,-0.5}};

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phiC1(int p, double x, double* dphi)
{
	double	f = 0.0;
	double	df = 0.0;
	double	signx = 0.0;
	short	piece = 0;
	short	i = 0;

	if (x > 0.0)
		signx = 1.0;
	else if (x < 0.0)
		signx = -1.0;

	x = fabs(x);
	if (x < p/2)
	{
		piece = (short)floor(x);
		x = x - (piece + 1.0);
		f = c[piece][p-1];
		df = c[piece][p-1]*(p-1);
		for (i = p-2; i > 0; i--)
		{
			f = f*x + c[piece][i];
			df = df*x + c[piece][i]*(i);
		}
		f = f*x + c[piece][0];
		df = df*signx;
	}

	if (dphi != NULL)
		*dphi = df;
	return f;
}

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of interpolant
g2fg is (p + 1)-vector
*/
void compute_g2fgC1(short p, double* g2fg)
{
	short		n = 0;

	assert(g2fg != NULL);

	for (n = 0; n <= p; n++)
	{
		g2fg[n] = phiC1(p,n/2.0,NULL);
	}
}

//	End of file