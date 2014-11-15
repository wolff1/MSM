//-------|---------|---------|---------|---------|---------|---------|---------|
/*
b_spline.c - Routines for initializing, using, and uninitializing B-splines
*/

#include "b_spline.h"

//	EXTERNAL Methods
void b_spline_initialize(void* Interpolant)
{
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;

	assert(Bs != NULL);
	printf("\tB_SPLINE initialization!\n");

	//	Set up COMMON function pointers
	Bs->cmn.evaluate = &b_spline_evaluate;
	Bs->cmn.uninitialize = &b_spline_uninitialize;

	//	Set up the B_SPLINE interpolant
	b_spline_compute_omega(Bs);
	b_spline_compute_omega_prime(Bs);
	b_spline_compute_g2p(Bs);
	b_spline_compute_g2fg(Bs);
	b_spline_compute_g2g(Bs);
	b_spline_compute_tg2g(Bs);
}

void b_spline_evaluate(void* Interpolant)
{
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;

	assert(Bs != NULL);

}

void b_spline_uninitialize(void* Interpolant)
{
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;

	assert(Bs != NULL);
	printf("\tUn-initializing B_SPLINE!\n");

	dynfree(Bs->omega);
	dynfree(Bs->omegap);
	dynfree(Bs->cmn.g2p);
	dynfree(Bs->cmn.g2fg);
	dynfree(Bs->cmn.g2g);
	dynfree(Bs->cmn.tg2g);
	dynfree(Bs);
}

//	INTERNAL Methods
void b_spline_compute_omega(B_SPLINE* Bs)
{
}

void b_spline_compute_omega_prime(B_SPLINE* Bs)
{
}

void b_spline_compute_g2p(B_SPLINE* Bs)
{
}

void b_spline_compute_g2fg(B_SPLINE* Bs)
{
}

void b_spline_compute_g2g(B_SPLINE* Bs)
{
}

void b_spline_compute_tg2g(B_SPLINE* Bs)
{
}

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

//	INTERNAL Methods
/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(short p, double x, double* dphi)
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

//  Horner's rule version of phi
void new_phi(short p, double** g2p, short n, double* x, double* phi, double* dphi)
{
    short       i = 0;
    short       j = 0;
    short       p_2 = p/2;
    short       k = 0;
    double      xp = 0.0;
    double      s = 0.0;

    assert(phi != NULL);
    assert(dphi != NULL);

    for (i = 0; i < n; i++)
    {
        phi[i] = 0.0;
        dphi[i] = 0.0;
        s = (x[i] > 0.0 ? 1.0 : -1.0);
        xp = fabs(x[i]);

        if (xp < p_2)
        {
            k = floor(xp);
            xp = xp - k - 1;

            phi[i] = g2p[k][p-1];
            dphi[i] = g2p[k][p-1]*(p-1);
            for (j = p-2; j > 0; j--)
            {
                phi[i] = phi[i]*xp + g2p[k][j];
                dphi[i] = dphi[i]*xp + g2p[k][j]*j;
            }
            phi[i] = phi[i]*xp + g2p[k][0];
            dphi[i] = dphi[i]*s;
        }
    }
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
		num *= (p-n);	// [(p)(p-1)(p-2)...(p/2+1)]
		den *= n;		// [(p/2)(p/2-1)...(2)(1)]*[2^(p-1)]
	}
	g2fg[0] = (double) num / (double) den;

	// Compute [p choose p/2 + n] with 1 <= n <= p/2
	for (n = 1; n <= p_2; n++)
	{
		g2fg[n] = (g2fg[n-1]*(double)(p_2-n+1)) / (double)(p_2+n);
	}
}

// End of file