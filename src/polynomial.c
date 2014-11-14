//-------|---------|---------|---------|---------|---------|---------|---------|
/*
polynomial.c - 
*/

#include "polynomial.h"

/*
a is input array of length a_degree representing polynomial coefficients
b is input array of length b_degree representing polynomial coefficients
c is output array of length a_degree+b_degree+1 representing polynomial
	coefficients s.t. c = a*b
NOTE:	It is assumed that the polynomials are in the same base
		By convention, use a for the higher degree, b for lower degree
*/
void mpoly(short a_degree, double* a, short b_degree, double* b, double* c)
{
	short		i = 0;
	short		j = 0;

	assert(a_degree > 0);
	assert(a != NULL);
	assert(b_degree > 0);
	assert(b != NULL);
	assert(c != NULL);

	// Ensure memory for c polynomial is clean before we use += operator
	for (i = 0; i <= a_degree+b_degree; i++)
	{
		c[i] = 0.0;
	}

	// Multiply a times b and set product to c
	for (i = 0; i <= a_degree; i++)	// could be carefully parallelized
	{
		for (j = 0; j <= b_degree; j++)
		{
			c[i+j] += a[i]*b[j];
		}
	}
}

/*
multiply polynomials where result has restricted degree

Same as above except that c is not to exceed c_degree.
NOTE: It is assumed that c has length c_degree+1
*/
void mpolyr(short a_degree, double* a, short b_degree, double* b,
			short c_degree, double* c)
{
	short		i = 0;
	short		j = 0;

	assert(a_degree > 0);
	assert(a != NULL);
	assert(b_degree > 0);
	assert(b != NULL);
	assert(c_degree > 0 && c_degree <= a_degree+b_degree+1);
	assert(c != NULL);

	// Ensure memory for c polynomial is clean before we use += operator
	for (i = 0; i <= c_degree; i++)
	{
		c[i] = 0.0;
	}

	// Multiply a times b and set product to c
	for (i = 0; i <= a_degree; i++)	// could be carefully parallelized
	{
		for (j = 0; j <= MIN(c_degree-i,b_degree); j++)
		{
			c[i+j] += a[i]*b[j];
		}
	}
}

//	End of file