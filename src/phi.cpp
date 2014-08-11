//-------|---------|---------|---------|---------|---------|---------|---------|
/*
centered b-spline phi
*/

#include "stdafx.h"

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(int p, double x, double* dphi)
{
	double f = 0.0;
	double df = 0.0;
	double** Q = NULL;
	double** Qp = NULL;
	double u = 0.0;
	short k = 0;
	short v = 0;

	if (fabs(x) < (double)(p/2.0))
	{
		// Dynamically create memory for Q and Qp 2D arrays
		Q = dynarr_d(p,p);		// Q(u)
		Qp = dynarr_d(p,p);		// Q'(u)
		u = x + (double)p/2.0;	// b/c centered B-spline

		// Start with k = 1 (Q_1 is the indicator function on [0,1[)
		k = 1;
		for (v = 0; v <= p-k; v++)
		{
			Q[0][v] = ((u-(double)v) < 1.0 && (u-(double)v) >= 0.0);
		}

		// use recurrance to build up to k=p-1
		for (k = 2; k <= p; k++)
		{
			for (v = 0; v <= p-k; v++)
			{
				Qp[k-1][v] = Q[k-2][v] - Q[k-2][v+1];
				Q [k-1][v] = (k*Q[k-2][v+1] + (u-(double)v)*Qp[k-1][v])/(k-1);
			}
		}

		// Set return variables
		f = Q[p-1][0];
		df = Qp[p-1][0];

		// Free dynamically allocated memory for Q and Qp 2D arrays
		dynfree(Q[0]);		// Free the data buffer memory
		dynfree(Qp[0]);
		dynfree(Q);			// Free the pointer-pointer memory
		dynfree(Qp);
	}

	if (dphi != NULL)
		*dphi = df;
	return f;
}

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of B-spline
g2fg is (p/2 + 1)-vector
*/
void compute_g2fg(short p, double* g2fg)
{
	short		n = 0;
	short		p_2 = (short)p/2;
	long		num = p;			// numerator
	long		den = p_2;			// denominator

	assert(g2fg != NULL);

	// Account for 2^(1-p) in denominator once
	for (n = 1; n < p; n++)
	{
		den *= 2;
	}

	// Build initial state of "p choose p/2"
	for (n = 1; n < p_2; n++) // exclude p_2 b/c it is in num and den
	{
		num *= (p-n);	// [(p)(p-1)(p-2)...(p-p/2)]
		den *= n;		// [(p/2)(p/2-1)...(2)(1)]*[2^(p-1)]
	}
	g2fg[0] = (double) num / (double) den;

	// Compute [p choose p/2 + n] with 1 <= n <= p/2
	for (n = 1; n <= p_2; n++)
	{
		g2fg[n] = (g2fg[n-1]*(double)(p_2-n+1)) / (double)(p_2+n);
	}
}

/*** DRIVER FUNCTIONS BELOW ***/

/*
Tests phi and phi' without regard to method or scaling.
*/
void phi_test_all(void)
{
	int			i = 0;
	int			samples = 0;
	int			p = 0;
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;

	/**************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%d", &p);
	assert(p % 2 == 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));

	/**************************************************************************/
	// Test centered B-spline where p must be even
	for (i = 0; i <= samples; i++)
	{
		X[i] = (double)(-p/2.0) + ((double)i/(double)samples)*((double)p);
		F[i] = phi(p, X[i], &DF[i]);
/*
		printf("%02d:\tphi(%f) = %f\tphi'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Phi(x)
	printf("Plotting Phi(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, F);
	//	Show Phi'(x)
	printf("Plotting Phi'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	/**************************************************************************/
	// Free allocated memory
	dynfree(X);
	dynfree(F);
	dynfree(DF);
}

/*
Simple interface to output nesting coefficients to user.
*/
void print_nesting_coefficients(void)
{
	short		p = 0;
	double*		J = NULL;

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%hd", &p);
	assert(p % 2 == 0);

	J = (double*) dynvec(p/2+1,sizeof(double));
	assert(J != NULL);

	// Calculate the coefficients to display
	compute_g2fg(p, J);

	for (short i = 0; i <= p/2; i++)
	{
		printf("J[%hd] = %lf\n", i, J[i]);
	}

	// Free dynamically allocated memory
	dynfree(J);
}

//	End of file