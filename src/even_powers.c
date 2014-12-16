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

	//	Initialize COMMON members
	Ep->cmn.Size = sizeof(EVEN_POWERS);

	//	Initialize COMMON function pointers
	Ep->cmn.copy = &even_powers_copy;
	Ep->cmn.soften = &even_powers_soften;
	Ep->cmn.split = &even_powers_split;
	Ep->cmn.uninitialize = &even_powers_uninitialize;

	// Initialize the EVEN_POWERS softener
	even_powers_compute_p2p(Ep);
}

void even_powers_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	Copy softening function coefficients
	((EVEN_POWERS*)Dst)->cmn.p2p = (double*) dynvec(((EVEN_POWERS*)Dst)->cmn.k+1, sizeof(double));
	memcpy(((EVEN_POWERS*)Dst)->cmn.p2p, ((EVEN_POWERS*)Src)->cmn.p2p, (((EVEN_POWERS*)Dst)->cmn.k+1)*sizeof(double));
}

void even_powers_soften(void* Softener, long Len, double* X, double* F, double* DF)
{	// aka "gamma"
	double			XX = 0.0;
	double*			c = NULL;
	long			i = 0;
	long			j = 0;
	short			k = 0;
	EVEN_POWERS*	Ep = (EVEN_POWERS*) Softener;

	assert(Ep != NULL);
	if (Len > 0)
	{
		assert(X != NULL);
		assert(F != NULL);
		assert(DF != NULL);
	}

	k = Ep->cmn.k;
	c = Ep->cmn.p2p;

	for (i = 0; i < Len; i++)
	{
		// gamma is a spline which becomes f(x) = 1/x for x >= 1.0
		if (X[i] < 1.0)
		{
			XX = X[i]*X[i];
			// Use Horner's rule to evaluate polynomial
			F[i] = c[k];			// Even powers
			DF[i] = c[k]*(2*k);	// Odd powers
			for (j = k-1; j >= 1; j--)
			{
				F[i] = F[i]*XX + c[j];
				DF[i] = DF[i]*XX + c[j]*(2*j);
			}
			F[i] = F[i]*XX + c[0];
			DF[i] = DF[i]*X[i];
		}
		else
		{
			F[i] = 1.0/X[i];
			DF[i] = -F[i]*F[i];	//	-f*f = -1.0/(x*x);
		}
	}
}

void even_powers_split(void* Softener, long Len, double* X, double* F, double* DF)
{	//	aka "theta"
	long			i = 0;
	double*			X_2 = NULL;
	double*			F_2 = NULL;
	double*			DF_2 = NULL;

	assert(Softener != NULL);
	assert(X != NULL);
	assert(F != NULL);
	assert(DF != NULL);

	X_2 = (double*) dynvec(Len, sizeof(double));
	F_2 = (double*) dynvec(Len, sizeof(double));
	DF_2 = (double*) dynvec(Len, sizeof(double));

	for (i = 0; i < Len; i++)
	{
		X_2[i] = 0.5*X[i];
	}

	even_powers_soften(Softener, Len, X, F, DF);
	even_powers_soften(Softener, Len, X_2, F_2, DF_2);

	for (i = 0; i < Len; i++)
	{
		F[i] = F[i] - 0.5*F_2[i];
		DF[i] = DF[i] - 0.25*DF_2[i];
	}

	dynfree(X_2);
	dynfree(F_2);
	dynfree(DF_2);
}

void even_powers_uninitialize(void* Softener)
{
	EVEN_POWERS*		Ep = (EVEN_POWERS*) Softener;
	assert(Ep != NULL);
	printf("\tUn-initializing EVEN_POWERS!\n");

	//	Free dynamically allocated memory for softening function coefficients
	dynfree(Ep->cmn.p2p);
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

// End of file