//-------|---------|---------|---------|---------|---------|---------|---------|
/*
b_spline.c - Routines related to operator(s)
*/

#include "b_spline.h"

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
Convert an operator from delta^2 to shifts: E-2+E

cd -> "central difference" operator (input)
s  -> "shifts" operator (output)

NOTE:
	- Both arrays should be length cd_degree+1
	- This routine does not initialize s
		- Instead it operates as s = s + ...

The kth term of (E-2+E)^n is given by:
s_k = ((-1)^(mod(n+k,2)))*nchoosek(2*n,k+n)
*/
void convert_delta2_to_shifts(short cd_degree, const double* cd, double* s)
{
	short		n = 0;
	short		k = 0;
	double		f = 1.0;	//	(2n)!
	double		nf = 1.0;	//	(n)!
	double		g = 1.0;	// coefficient for a term due to conversion
	double		sign1 = 1.0;
	double		sign2 = 1.0;

	assert(cd_degree > -1);
	assert(cd != NULL);
	assert(s != NULL);

	//	Convert cd(\delta^2) to s(E)
	for (n = 0; n <= cd_degree; n++)
	{
		g = f/(nf*nf);
		sign2 = sign1;
		s[0] = s[0] + sign2*cd[n]*g;
		for (k = 1; k <= n; k++)
		{
			sign2 = -sign2;
			g = ((n-k+1.0)*g) / ((double)(n+k));
			s[k] = s[k] + sign2*cd[n]*g;
		}

		// Update the (2n)! part of the combination for the next iteration
		f = (2.0*n+2.0)*(2.0*n+1.0)*f;
		nf = (n+1.0)*nf;
		sign1 = -sign1;
	}
}

/*
compute omega values
NOTE:
	B has degree p-1
*/
void compute_omega(short p, short n, double* omega)
{
	short		i = 0;
	short		p_2 = p/2;
	double*		B = NULL;
	double*		Psi = NULL;

	assert(p%2 == 0);
	assert(n >= 0);
	assert(omega != NULL);

	// Dynamically allocate memory
	B = (double*) dynvec(p_2, sizeof(double));
	Psi = (double*) dynvec(1,sizeof(double));

	//	Compute B(E)
	for (i = 0; i < p_2; i++)
	{
		B[i] = phi(p,(double)i,NULL);
	}

	//	Set up Psi
	Psi[0] = 1.0;

	//	Solve for omega
    bibst_lss(-1, pow(2.0,-53), p_2, B, 1, 1, Psi, n, omega);
/*
    for (i = 0; i < p_2; i++)
    {
        printf("c[%d] = %f\n", i, omega[i]);
    }
*/
	// Free dynamically allocated memory
	dynfree(B);
	dynfree(Psi);
}

/*
compute omega' values
NOTE:
	B has degree p-1
	B2 has degree 2p-2
	A2 has degree p-1
	omega' has degree p/2+mu
*/
void compute_omega_prime(short p, short mu, double* omegap)
{
	short		i = 0;
	short		p_2 = p/2;
	double*		B = NULL;
	double*		B2 = NULL;
    double*     B2E = NULL;
	double*		A2 = NULL; // A^2 ~= B^{-2}
	double*		Ctmp = NULL;
	double*		C = NULL;
	double*		CE = NULL;
	double*		AE = NULL;
    double*     c = NULL;
	double*		pd = NULL;
	double*		pE = NULL;
	double*		c_full = NULL;
	double*		p_full = NULL;
	double*		high_terms = NULL;

	assert(p%2 == 0);
	assert(mu >= 0);
	assert(omegap != NULL);

	// Dynamically allocate memory
	B = (double*) dynvec(p_2, sizeof(double));
	B2 = (double*) dynvec(p-1, sizeof(double));
    B2E = (double*) dynvec(p-1, sizeof(double));
	A2 = (double*) dynvec(p_2, sizeof(double));
	Ctmp = (double*) dynvec(p+p_2-2, sizeof(double));
	C = (double*) dynvec(p-2, sizeof(double));
	CE = (double*) dynvec(p-2, sizeof(double));
    c = (double*) dynvec(mu+1, sizeof(double));
	pd = (double*) dynvec(p_2+1,sizeof(double));
	pE = (double*) dynvec(p_2+1,sizeof(double));
	c_full = (double*) dynvec(2*mu+1,sizeof(double));
	p_full = (double*) dynvec(p+1,sizeof(double));
	high_terms = (double*) dynvec(p+2*mu+1,sizeof(double));

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
	//	compute (A^2)(B^2) in order to get -C
	mpoly(p-2, B2, p_2-1, A2, Ctmp);
	//	keep only (delta^2)^{p/2} and higher
	for (i = p_2; i <= p+p_2-3; i++)
	{
		C[i-p_2] = -Ctmp[i];
	}

	//	convert C(delta^2) to C(E)
	convert_delta2_to_shifts(p-3, C, CE);
/*
	for (i = 0; i <= p-3; i++)
	{
		printf("CE[%d] = %f\n", i, CE[i]);
	}
	printf("\n");
*/
	//	convert A^2(delta^2) to A(E) = omega'
	convert_delta2_to_shifts(p_2-1, A2, omegap);
/*
	for (i = 0; i <= p_2-1; i++)
	{
		printf("omega'[%d] = %f\n", i, omegap[i]);
	}
	printf("\n");
*/
	// Solve bi-infinite linear system involving B^2(E) and C(E)
    convert_delta2_to_shifts(p-2, B2, B2E);
    bibst_lss(-1, pow(2.0,-53), p-1, B2E, p-2, p-2, CE, mu+1, c);
/*
    for (i = 0; i <= p+mu; i++)
    {
        printf("c[%d] = %f\n", i, c[i]);
    }
*/
	// use solution, c, to find (\delta^2)^{p/2} \sum c_m E^m in E basis
	pd[p_2] = 1.0;
	convert_delta2_to_shifts(p_2,pd,pE);
/*
	display_vector_d(pE,p_2+1);
*/
	// convolve (c,pE) after forming full, symmetric operator
	c_full[mu] = c[0];
	for (i = 1; i <= mu; i++)
	{
		c_full[mu+i] = c[i];
		c_full[mu-i] = c[i];
	}
//display_vector_d(c_full,2*mu+1);
	p_full[p_2] = pE[0];
	for (i = 1; i <= p_2; i++)
	{
		p_full[p_2+i] = pE[i];
		p_full[p_2-i] = pE[i];
	}
//display_vector_d(p_full,p+1);
	mpoly(2*mu, c_full, p, p_full, high_terms);
//display_vector_d(high_terms,p+2*mu+1);

	//	add previous step to omega'
	for (i = 0; i <= p_2+mu; i++)
	{
		omegap[i] += high_terms[p_2+mu+i];
	}

	// Free dynamically allocated memory
	dynfree(B);
	dynfree(B2);
    dynfree(B2E);
	dynfree(A2);
	dynfree(Ctmp);
	dynfree(C);
	dynfree(CE);
    dynfree(c);
	dynfree(pd);
	dynfree(pE);
	dynfree(c_full);
	dynfree(p_full);
	dynfree(high_terms);
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
driver to test building the omega values
*/
void test_omega(void)
{
	short		i = 0;
	short		p = 0;
	short		p_2 = 0;
	short		n = 0;
	double*		omega = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);
	p_2 = p/2;
	printf("Please enter n: ");
	scanf("%hd", &n);

	// Dynamically allocate memory
	omega = (double*) dynvec(n, sizeof(double));

	compute_omega(p, n, omega);

	printf("\n");
	for (i = 0; i < n; i++)
	{
		printf("%02hd: %f\n", i, omega[i]);
	}
	printf("\n");

	// Free dynamically allocated memory
	dynfree(omega);
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

	printf("\n");
	for (i = 0; i <= p_2+mu; i++)
	{
		printf("%02hd: %f\n", i, omegap[i]);
	}
	printf("\n");

	// Free dynamically allocated memory
	dynfree(omegap);
}

/*
driver to test conversion from delta^2 to shifts
*/
void test_convert_to_shifts(void)
{
	double*		cd = NULL;
	double*		s = NULL;
	short		degree = 2;
	short		i = 0;
	short		j = 0;

	for (i = 0; i <= 4; i++)
	{
		// Allocate
		cd = (double*) dynvec(i+1, sizeof(double));
		s  = (double*) dynvec(i+1, sizeof(double));

		// set one term and display
		cd[i] = 1.0;
		for (j = 0; j <= i; j++)
		{
			printf("cd[%02hd] = %f\n", j, cd[j]);
		}
		printf("\n");

		// Convert cd to s
		convert_delta2_to_shifts(i, cd, s);

		// display result
		for (j = 0; j <= i; j++)
		{
			printf("s[%02hd] = %f\n", j, s[j]);
		}
		printf("\n----------------------------------\n");

		// Free
		dynfree(cd);
		dynfree(s);
	}
}

// End of file