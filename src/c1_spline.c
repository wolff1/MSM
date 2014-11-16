//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.c - C1 phi (CUBIC ONLY)
*/

#include "c1_spline.h"

//	EXTERNAL Methods
void c1_spline_initialize(void* Interpolant, MSM_PARAMETERS* MsmParams)
{
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;
	assert(C1 != NULL);
	assert(MsmParams != NULL);
	printf("\tC1_SPLINE initialization!\n");

	//	Set up COMMON function pointers
	C1->cmn.compute_g2g = &c1_spline_compute_g2g;
	C1->cmn.compute_tg2g = &c1_spline_compute_tg2g;
	C1->cmn.evaluate = &c1_spline_evaluate;
	C1->cmn.uninitialize = &c1_spline_uninitialize;

	//	Set up the C1_SPLINE interpolant
	c1_spline_compute_g2p(C1);
	c1_spline_compute_g2fg(C1);
}

void c1_spline_compute_g2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams)
{
	long			Size = 0;
	double			Scale = 1.0;
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;

	assert(C1 != NULL);
	assert(Softener != NULL);
	assert(MsmParams != NULL);
	printf("\tC1_SPLINE compute_g2g\n");

	//	Pre-processing (Intermediate levels)
	Size = (long) ceil(2.0*MsmParams->alpha);
	Scale = MsmParams->h/MsmParams->a;

	//	Compute K^l sequence (defined by theta function) for intermediate level grid(s)
	stencil_initialize(&C1->cmn.g2g, Size, STENCIL_SHAPE_SPHERE);
	stencil_populate(C1->cmn.g2g, Softener, STENCIL_FUNCTION_TYPE_THETA, Scale);
//	stencil_display(C1->cmn.g2g, MsmParams->h/MsmParams->a);
}

void c1_spline_compute_tg2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams)
{
	long			Size = 0;
	double			Scale = 1.0;
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;

	assert(C1 != NULL);
	assert(Softener != NULL);
	assert(MsmParams != NULL);
	printf("\tC1_SPLINE compute_tg2g\n");

	//	Pre-processing (Top level)
	Size = (long) ceil(MsmParams->D);
	Scale = MsmParams->h/MsmParams->a;

	//	Compute K^L sequence (defined by gamma function) for top level grid
	stencil_initialize(&C1->cmn.tg2g, Size, STENCIL_SHAPE_CUBE);
	stencil_populate(C1->cmn.tg2g, Softener, STENCIL_FUNCTION_TYPE_GAMMA, Scale);
//	stencil_display(C1->cmn.tg2g, MsmParams->h/MsmParams->a);
}

void c1_spline_evaluate(void* Interpolant, long Len, double* X, double* F, double* DF)
{
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;
	assert(C1 != NULL);
	printf("\tC1_SPLINE evaluate\n");
}

void c1_spline_uninitialize(void* Interpolant)
{
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;
	assert(C1 != NULL);
	printf("\tUn-initializing C1_SPLINE!\n");

	//	Free the dynamically allocated memory
	stencil_free(C1->cmn.g2g);
	stencil_free(C1->cmn.tg2g);
	//dynfree(C1->cmn.g2p[0]);
	//dynfree(C1->cmn.g2p);
	//dynfree(C1->cmn.g2fg);
	//dynfree(C1->cmn.g2g);
	//dynfree(C1->cmn.tg2g);
	dynfree(C1);
}

//	INTERNAL Methods
void c1_spline_compute_g2p(C1_SPLINE* C1)
{
	assert(C1 != NULL);
	printf("\tC1_SPLINE compute_g2p\n");
}

void c1_spline_compute_g2fg(C1_SPLINE* C1)
{
	assert(C1 != NULL);
	printf("\tC1_SPLINE compute_g2fg\n");
}

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phiC1(int p, double x, double* dphi)
{
	// These are for the CUBIC "shifted-powers" version
	double c[2][4] = {{0.0, -0.5, 2.0, 1.5},
					  {0.0, 0.0, -0.5,-0.5}};

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