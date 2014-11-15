//-------|---------|---------|---------|---------|---------|---------|---------|
/*
even_powers.c - the softening function used to split the kernel
*/

#include "even_powers.h"

//	EXTERNAL Methods
void even_powers_initialize(void* Softener)
{
	EVEN_POWERS*		Ep = (EVEN_POWERS*) Softener;
	assert(Ep != NULL);
	printf("\tEVEN_POWERS initialization!\n");

	//	Initialize COMMON function pointers
	Ep->cmn.soften = &even_powers_soften;
	Ep->cmn.split = &even_powers_split;
	Ep->cmn.uninitialize = &even_powers_uninitialize;

	// Initialize the EVEN_POWERS softener
	even_powers_compute_p2p(Ep);
}

void even_powers_soften(void* Softener, long Len, double* X, double* FX, double* DFX)
{	// aka "gamma"
	double			XX = 0.0;
	double*			c = NULL;
	long			i = 0;
	long			j = 0;
	short			k = 0;
	EVEN_POWERS*	Ep = (EVEN_POWERS*) Softener;

	assert(Ep != NULL);
	assert(X != NULL);
	assert(FX != NULL);
	assert(DFX != NULL);

	k = Ep->cmn.k;
	c = Ep->cmn.p2p;

	for (i = 0; i < Len; i++)
	{
		// gamma is a spline which becomes f(x) = 1/x for x >= 1.0
		if (X[i] >= 1.0)
		{
			FX[i] = 1.0/X[i];
			DFX[i] = -FX[i]*FX[i];	//	-f*f = -1.0/(x*x);
		}
		else
		{
			XX = X[i]*X[i];
			// Use Horner's rule to evaluate polynomial
			FX[i] = c[k];			// Even powers
			DFX[i] = c[k]*(2*k);	// Odd powers
			for (j = k-1; j >= 1; j--)
			{
				FX[i] = FX[i]*XX + c[j];
				DFX[i] = DFX[i]*XX + c[j]*(2*j);
			}
			FX[i] = FX[i]*XX + c[0];
			DFX[i] = DFX[i]*X[i];
		}
	}
}

void even_powers_split(void* Softener, long Len, double* X, double* FX, double* DFX)
{	//	aka "theta"
	long			i = 0;
	double*			X_2 = NULL;
	double*			FX_2 = NULL;
	double*			DFX_2 = NULL;

	assert(Softener != NULL);
	assert(X != NULL);
	assert(FX != NULL);
	assert(DFX != NULL);

	X_2 = (double*) dynvec(Len, sizeof(double));
	FX_2 = (double*) dynvec(Len, sizeof(double));
	DFX_2 = (double*) dynvec(Len, sizeof(double));

	for (i = 0; i < Len; i++)
	{
		X_2[i] = 0.5*X[i];
	}

	even_powers_soften(Softener, Len, X, FX, DFX);
	even_powers_soften(Softener, Len, X_2, FX_2, DFX_2);

	for (i = 0; i < Len; i++)
	{
		FX[i] = FX[i] - 0.5*FX_2[i];
		DFX[i] = DFX[i] - 0.25*DFX_2[i];
	}

	dynfree(X_2);
	dynfree(FX_2);
	dynfree(DFX_2);
}

void even_powers_uninitialize(void* Softener)
{
	EVEN_POWERS*		Ep = (EVEN_POWERS*) Softener;
	assert(Ep != NULL);
	printf("\tUn-initializing EVEN_POWERS!\n");

	//	Free dynamically allocated memory for softening function coefficients
	dynfree(Ep->cmn.p2p);

	//	Free dynamically allocated memory for EVEN_POWERS object
	dynfree(Ep);
}

//	INTERNAL Methods
void even_powers_compute_p2p(EVEN_POWERS* Ep)
{
	short		k = Ep->cmn.k;
	double**	A = NULL;	// Coefficient matrix
	double*		b = NULL;	// rhs vector first, then solution vector
	int			i = 0;
	int			j = 0;
	lapack_int	rc = 0;
	lapack_int*	piv = NULL;

	// Allocate memory for A and b
	A = dynarr_d(k+1, k+1);
	b = (double*) dynvec(k+1,sizeof(double));

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

	// Solve system Ax = b
	// NOTE: b is overwritten with x, A is overwritten with LU
	piv = (lapack_int*) dynvec(k+1,sizeof(lapack_int));
	rc = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) k+1, (lapack_int) 1,
						A[0], (lapack_int) k+1, piv, b, (lapack_int) 1);
	dynfree(piv);
	assert(rc == 0); // zero is SUCCESS

	//	b is now solution vector (i.e. the coefficients of the softening function)
	Ep->cmn.p2p = b;

	// Free allocated memory
	dynfree(A[0]);
	dynfree(A);
}

//*****************************************************//
//*****DELETE EVERYTHING FROM HERE DOWN EVENTUALLY*****//
//*****************************************************//

/*
Solves for coefficients of gamma
k is degree of continuity
b will return k+1 vector of coefficients
*/
void gamma_init(short k, double* x)
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
double _gamma(double *c, short k, double x, double* dgamma)
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
*/
double theta(double *c, short k, double x, double* dtheta)
{
	double f = 0.0;
	double df = 0.0;
	double df2 = 0.0;

    f = _gamma(c, k, x, &df) - 0.5*_gamma(c, k, 0.5*x, &df2);
    df = df - 0.25*df2;

	if (dtheta != NULL)
		*dtheta = df;
	return f;
}

/*
Evaluate theta* and theta*' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta_star(double *c, short k, double x, double* dtheta_star)
{
    double f = 0.0;
    double df = 0.0;

    f = 1.0/x - _gamma(c, k, x, &df);
    df = -1.0/(x*x) - df;

    if (dtheta_star != NULL)
        *dtheta_star = df;
    return f;
}

/*
Evaluate theta and theta' at position x as single polynomial
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double thetap(double *c, short k, double x, double* dtheta)
{
	double	f = 0.0;
	double	df = 0.0;
	double	exp2 = 0.0;
	double	xx = x*x;
	short	i = 0;

	if (x <= 1.0)
	{
		exp2 = pow(2.0,2*k+1);
		// Use Horner's rule to evaluate polynomial
		f = ((exp2-1.0)/exp2)*c[k];			// Even powers
		df = ((exp2-1.0)/exp2)*c[k]*(2*k);	// Odd powers
		for (i = k-1; i >= 1; i--)
		{
			exp2 *= 0.25;//exp2 = pow(2.0,2*i+1);
			f = f*xx + ((exp2-1.0)/exp2)*c[i];
			df = df*xx + ((exp2-1.0)/exp2)*c[i]*(2*i);
		}
		f = f*xx + (0.5)*c[0];
		df = df*x;
	}
	else if (x < 2.0)
	{
		f = 1.0/x - 0.5*_gamma(c,k,0.5*x,&df);
		df = -1.0/xx - 0.25*df;
	}

	if (dtheta != NULL)
		*dtheta = df;
	return f;
}

// End of file