//-------|---------|---------|---------|---------|---------|---------|---------|
/*
Routines related to operator(s)
*/
#include "stdafx.h"

/*
void compute_blurring_operator(short degree)
where
	degree is the max degree of the operator in s^2 and
	B is the blurring operator.

To calculate blurring operator, B_p:
c_n = 1/a_0 (b_n - [a_1*c_{n-1} + a_2*c_{n-2} + ... + a_n*c_0])

NOTE:
	* we do not need to store b at all
	* factorials are not calculated anew in every iteration of the loops
		* f keeps a running product that is equivalent to factorial()
	* inner loop only touches columns which will change for given n,k
		* inner loop could be parallelized, although probably not worth it
	* for more details, go to
		http://pmsm.cs.purdue.edu/tiki-index.php?page=omegaprime
		or
		http://en.wikipedia.org/wiki/Formal_power_series#Dividing_series
*/
void compute_blurring_operator(short degree, double* B)
{
	double*		a = NULL;
	double**	c = NULL;
	double		f = 1.0;	// factorial
	short		n = 0;
	short		k = 0;
	short		i = 0;

	assert(B != NULL);

	//	Allocate memory dynamically for computation
	a = (double*) dynvec(degree+1,sizeof(double));
	c = (double**) dynarr_d(degree,degree+1);

	// All but last iteration
	a[0] = 1.0;
	c[0][0] = 1.0;			// b_0 / a_0
	for (n = 1; n < degree; n++)
	{
		a[n] = -1.0 / (n*f);
		f = (2.0*n+1.0)*(2.0*n)*f;	// Keep running factorial
		c[n][n] = 1.0 / f;			// b_n
		for (k = 1; k <= n; k++)
		{
			for (i = 0; i <= n-k; i++)
			{
				c[n][i+k-1] = c[n][i+k-1] - a[k]*c[n-k][i];
			}
		}
	}

	// Last iteration, use B instead of c
	n = degree;
	a[n] = -1.0 / (n*f);
	f = (2.0*n+1.0)*(2.0*n)*f;	// Keep running factorial
	B[n] = 1.0 / f;				// c_n = b_n
	for (k = 1; k <= n; k++)
	{
		for (i = 0; i <= n-k; i++)
		{
			B[i+k-1] = B[i+k-1] - a[k]*c[n-k][i];
		}
	}

	//	Free dynamically allocated memory
	dynfree(a);
	dynfree(c[0]);
	dynfree(c);
}

/*
driver to test compute_blurring_operator
*/
void test_blurring_operator(void)
{
	short	p = 0;
	short	k = 0;
	double*	B = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);

	// Dynamically allocate memory for blurring operator
	B = (double*) dynvec(p,sizeof(double));

	compute_blurring_operator(p-1, B);

	// Print c_n
	printf("B = (%f", B[0]);
	for (k = 1; k < p; k++)
	{
		printf(", %f", B[k]);
	}
	printf(")\n");

	// Free dynamically allocated memory
	dynfree(B);
}
