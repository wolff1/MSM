//-------|---------|---------|---------|---------|---------|---------|---------|
/*
gamma.cpp - the softening function used to split the kernel
*/

#include "stdafx.h"

/*
Solves for coefficients of gamma
k is degree of continuity
b will return k+1 vector of coefficients
*/
void gamma_init(int k, double* x)
{
	double** A = NULL;	// Coefficient matrix
	double* b = x;		// Alias output parameter x
	int i = 0;
	int j = 0;
	lapack_int rc = 0;
	lapack_int* piv = NULL;

	// b will first act as RHS, then return coefficients to caller
	assert(b != NULL);

	// Allocate memory for A
	A = dynarr_d(k+1, k+1);

	// Build row zero of A and b
	b[0] = 1.0;
	for (i = 0; i <= k; i++)
	{
		A[0][i] = 1.0;
	}

	// Build rows 1 through k of A and b
	for (i = 1; i <= k; i++)
	{
		b[i] = -i*b[i-1];
		// cols (starting with diagonal or subdiagonal)
		for (j = ceil(i/2.0); j <= k; j++)
		{
			// NOTE: 2j-i+1 is the power of the term of one less derivative
			A[i][j] = (2*j-i+1)*A[i-1][j];
		}
	}
/*
	// Display A
	printf("A = \n");
	display_dynarr_d(A,k+1, k+1);

	// Display b
	printf("b = \n");
	display_vector_d(b, k+1);
*/
	// Solve system Ax = b
	// NOTE: b is overwritten with x, A is overwritten with LU
	piv = (lapack_int*) dynvec(k+1,sizeof(lapack_int));
	rc = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) k+1, (lapack_int) 1,
						A[0], (lapack_int) k+1, piv, b, (lapack_int) 1);
	dynfree(piv);
	assert(rc == 0); // zero is SUCCESS
/*
	// Display c
	printf("c = \n");
	display_vector_d(b, k+1);
*/
	// Free allocated memory
	dynfree(A[0]);
	dynfree(A);
}

/*
Evaluate gamma and gamma' at (positive) position x
c is coefficient vector for gamma
*/
double gamma(double *c, int k, double x, double* dgamma)
{
	double f = 0.0;
	double df = 0.0;
	double xx = x*x;
	int i = 0;

	// gamma is a spline which becomes f(x) = 1/x for x >= 1.0
	if (x >= 1.0)
	{
		f = 1.0/x;
		df = -f*f;	//	-f*f = -1.0/(x*x);
	}
	else
	{
		// Use Horner's rule to evaluate polynomial
		f = c[k];			// Even powers
		df = c[k]*(2*k);	// Odd powers
		for (i = k-1; i >= 1; i--)
		{
			f = f*xx + c[i];
			df = df*xx + c[i]*(2*i);
		}
		f = f*xx + c[0];
		df = df*x;
	}

	if (dgamma != NULL)
		*dgamma = df;
	return f;
}

/*
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
long_range is flag to designate short (theta*) or long (theta) range 
*/
double theta(double *c, int k, double x, double* dtheta, int long_range)
{
	double f = 0.0;
	double df = 0.0;
	double f2 = 0.0;
	double df2 = 0.0;

/*
	f = gamma(c, k, x, &df);

	if (long_range == 0)
	{
		//f = 1.0/x - f;
		//df = -1.0/(x*x) - df;
		f = -f + 1.0/x;
		df = -df - 1.0/(x*x);
	}
	else
	{
		f2 = 0.5*gamma(c, k, 0.5*x, &df2);
		df2 = 0.25*df2;

		f = f - f2;
		df = df - df2;
	}
*/

	// Below should be the same as above, but without conditionals
	f = gamma(c, k, x, &df);
	f2 = 0.5*gamma(c, k, 0.5*x, &df2);

	// Compute theta(x)
	f = (double)(long_range)*(f - f2) +					// theta
		(double)(1 - long_range)*(-f + 1.0/x);			// theta*

	// Compute theta'(x)
	df = (double)(long_range)*(df - 0.25*df2) +			// theta
		 (double)(1 - long_range)*(-df - 1.0/(x*x));	// theta*

	if (dtheta != NULL)
		*dtheta = df;
	return f;
}

// End of file