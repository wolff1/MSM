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

	assert(degree > 0);
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
Use O^{-1} = 1 / O = 1 / 1 - D -> D = 1 - O
Then use geometric series
O^{-1} = 1 + D + D^2 + D^3 + ... + D^{I_degree}

NOTE:	O ("oh", not zero) is the input operator
		I ("eye") is the inverse, I = O^{-1}
*/
void compute_operator_inverse(short O_degree, const double* O,
							  short I_degree, double* I)
{
	short		i = 0;
	short		imax = MIN(O_degree,I_degree);
	short		k = 0;
	double*		D = NULL;
	double*		ID = NULL;

	assert(O_degree > 0);
	assert(O != NULL);
	assert(I_degree > 0);
	assert(I != NULL);

	// Dynamically allocate memory
	D = (double*) dynvec(O_degree+1, sizeof(double));
	ID = (double*) dynvec(I_degree+1, sizeof(double));

	// Initialize D = 1 - O
	D[0] = 0.0;	// 1 - O[0] where O[0] is 1
	for (i = 1; i <= O_degree; i++)
	{
		D[i] = -O[i];
	}

	// Ensure that memory for I is clean
	for (i = 0; i <= I_degree; i++)
	{
		I[i] = 0.0;
	}

	//	Set I using Horner's Rule
	//	First, I = 1 + D
	I[0] = 1.0;
	for (i = 1; i <= imax; i++)
	{
		I[i] = D[i];
	}

	// Then, iterate for each necessary power of D
	for (k = I_degree-1; k > 0; k--)
	{
		//I = I*D + 1;
		mpolyr(I_degree, I, O_degree, D, I_degree, ID);
		for (i = 0; i <= I_degree; i++)
		{
			I[i] = ID[i];
		}
		I[0] += 1.0;
	}

	// Free dynamically allocated memory
	dynfree(D);
	dynfree(ID);
}

/*
compute omega' values
NOTE:
	B has degree p-1
	B2 has degree 2p-2
	A2 has degree p-1
*/
void compute_omega_prime(short p, short mu, double* omegap)
{
	short		i = 0;
	short		p_2 = p/2;
	double*		B = NULL;
	double*		B2 = NULL;
	double*		A2 = NULL; // A^2 ~= B^{-2}

	assert(p%2 == 0);
	assert(mu >= 0);
	assert(omegap != NULL);

	// Dynamically allocate memory
	B = (double*) dynvec(p_2, sizeof(double));
	B2 = (double*) dynvec(p-1, sizeof(double));
	A2 = (double*) dynvec(p_2, sizeof(double));

	//	Compute B_p
	compute_blurring_operator(p_2-1, B);
/*
	printf("Blurring (1D):\n");
	for (i = 0; i < p_2; i++)
	{
		printf("%02hd - %f\n", i, B[i]);
	}
	printf("\n");
*/
	//	Compute B^2
	mpoly(p_2-1, B, p_2-1, B, B2);
/*
	printf("Blurring (2D):\n");
	for (i = 0; i < p-1; i++)
	{
		printf("%02hd - %f\n", i, B2[i]);
	}
	printf("\n");
*/
	//	Compute B^{-2} = 1 / B^2 = 1 / 1 - D => D = 1 - B^2
	compute_operator_inverse(p-2, B2, p_2-1, A2);
/*
	printf("Antiblurring (2D):\n");
	for (i = 0; i <= p_2-1; i++)
	{
		printf("%02hd - %f\n", i, A2[i]);
	}
	printf("\n");
*/
	// Free dynamically allocated memory
	dynfree(B);
	dynfree(B2);
	dynfree(A2);
}

/*
driver to test compute_blurring_operator
*/
void test_blurring_operator(void)
{
	short	p = 0;
	short	p_2 = 0;
	short	k = 0;
	double*	B = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);
	p_2 = p/2;

	// Dynamically allocate memory for blurring operator
	B = (double*) dynvec(p_2,sizeof(double));

	compute_blurring_operator(p_2-1, B);

	// Print c_n
	printf("B = (%f", B[0]);
	for (k = 1; k < p_2; k++)
	{
		printf(", %f", B[k]);
	}
	printf(")\n");

	// Free dynamically allocated memory
	dynfree(B);
}

/*
driver to test polynomial multiplication
*/
void test_mpoly(void)
{
	double*		a = NULL;
	double*		b = NULL;
	double*		c = NULL;
	short		i = 0;
	short		d1 = 0;
	short		d2 = 0;
	short		d3 = 0;

	//	Get degrees from user
	printf("Please enter degree 1: ");
	scanf("%hd", &d1);
	printf("Please enter degree 2: ");
	scanf("%hd", &d2);
	d3 = d1+d2;

	// Dynamically allocate memory
	a = (double*) dynvec(d1+1, sizeof(double));
	b = (double*) dynvec(d2+1, sizeof(double));
	c = (double*) dynvec(d3+1, sizeof(double));

	// Initialize a
	for (i = 0; i <= d1; i++)
	{
		a[i] = (double)i + 1.0;
	}

	// Initialize b
	for (i = 0; i <= d2; i++)
	{
		b[i] = (double)-i - 1.0;
	}

	mpoly(d1, a, d2, b, c);

	printf("\n");
	for (i = 0; i <= d3; i++)
	{
		printf("%f ", c[i]);
	}
	printf("\n");

	// Free dynamically allocated memory
	dynfree(a);
	dynfree(b);
	dynfree(c);
}

/*
driver to test building the omega' values
*/
void test_omegap(void)
{
	short		i = 0;
	short		p = 0;
	short		p_2 = 0;
	short		mu = 0;
	double*		omegap = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);
	p_2 = p/2;
	printf("Please enter mu: ");
	scanf("%hd", &mu);

	// Dynamically allocate memory
	omegap = (double*) dynvec(p_2+mu+1, sizeof(double));

	compute_omega_prime(p, mu, omegap);
/*
	printf("\n");
	for (i = 0; i < p_2+mu+1; i++)
	{
		printf("%02hd - %f\n", i, omegap[i]);
	}
	printf("\n");
*/
	// Free dynamically allocated memory
	dynfree(omegap);
}

// End of file