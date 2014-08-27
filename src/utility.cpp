//-------|---------|---------|---------|---------|---------|---------|---------|
/*
utility functions
*/

#include "stdafx.h"

/*
NOTE: See https://software.intel.com/en-us/node/433552 for data alignment
NOTE: See https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/
		GUID-637284D3-4D1F-4D6C-9509-382CB2DD1A3D.htm for routines involving
		linear systems
*/

/*
Allocate and zero memory for 2D array of type double
*/
double** dynarr_d(int rows, int cols)
{
	int i = 0;
	double** ptr = NULL;
	double* buf = NULL;

	// Allocate pointer buffer
	//ptr = (double**) calloc(rows, sizeof(double*));
	ptr = (double**) mkl_malloc(rows*sizeof(double*), MEM_ALIGN);
	assert(ptr != NULL);
	memset(ptr, 0, rows*sizeof(double*));

	// Allocate data buffer
	buf = (double*) dynvec(rows*cols, sizeof(double));

	// Assign ptr to point into buffer
	for (i = 0; i < rows; i++)
	{
		ptr[i] = buf + (i*cols);
	}

	return ptr;
}

/*
Allocate and zero memory for vector (1D array) whose elements
have size "size"
*/
void* dynvec(int rows, size_t size)
{
	void* ptr = NULL;

	// Allocate pointer buffer
	//ptr = (double*) calloc(rows, sizeof(double));
	ptr = mkl_malloc(rows*size, MEM_ALIGN);
	assert(ptr != NULL);
	memset(ptr, 0, rows*size);

	return ptr;
}

/*
Free allocated memory
*/
void dynfree(void* ptr)
{
	//free(ptr);	// Consider mkl_free()
	mkl_free(ptr);
}

/*
Display elements of 2D array
*/
void display_dynarr_d(double** A, int rows, int cols)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			printf("%f ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*
Display elements of vector
*/
void display_vector_d(double* x, int rows)
{
	int i = 0;

	for (i = 0; i < rows; i++)
	{
		printf("%f\n", x[i]);
	}
	printf("\n");
}

/*
a is input array of length a_degree representing polynomial coefficients
b is input array of length b_degree representing polynomial coefficients
c is output array of length a_degree+b_degree+1 representing polynomial coefficients
	s.t. c = a*b
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
void mpolyr(short a_degree, double* a, short b_degree, double* b, short c_degree, double* c)
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


/*
bibst -> bi-infinite, banded, symmetric, Toeplitz
lss -> linear system solver

A_len   -> bandwidth
A       -> band (diagonal to edge of band)
b_len   -> length of rhs vector b, gives # of equations to solve
b_nnz   -> the number of non-zero elements in b
b       -> rhs vector
x_len   -> length of the solution vector x
x       -> solution vector (Ax = b)
*/
void bibst_lss(long max_itr, double tol, short A_len, double* A, short b_len, short b_nnz, double* b, short x_len, double* x)
{
    short           i = 0;
    short           j = 0;
    short           k = 0;
    short           bw = A_len; //  Change A_len to bw in definition?
    short           imin = 0;   //  Used in infinite Cholesky
    double          ods = 0.0;  //  off-diagonal sum
    double          ds = 0.0;   //  on-diagonal sum
    double          norm = 0.0; //  Used in convergence check
    double**        G = NULL;   //  Lower triangular Cholesky factor
    double*         pr = NULL;  //  Previous Row (used in infinite routine)
    char            c[2] = {0,0};
    
    assert(A_len > 0);
    assert(A != NULL);
    assert(b_len > 0);
    assert(b_nnz > 0);
    assert(b != NULL);
    assert(x_len > 0);
    assert(x != NULL);

    //  Dynamically allocate memory for matrices and vectors
    G = (double**) dynarr_d(bw,bw);
    pr = (double*) dynvec(bw, sizeof(double));

    if (max_itr < 0)
        max_itr = 1000;

    //  Default values in G to be 1 on diagonal so initial divisions by
    //  diagonal element do not produce non-values
    for (i = 0; i < bw; i++)
    {
        G[i][i] = 1.0;
    }

/*
    Factor A = (G^T)(G) (Upper-Lower Cholesky)
        -> Use infinite Cholesky algorithm to get "converged" row
        -> Use finite Cholesky algorithm to get "unique" rows
*/

    for (j = 0; j < max_itr; j++)
    {
        ds = 0.0;
        imin = MAX(j - bw + 1, 0);
        for (i = imin; i < j; i++)
        {
            ods = 0.0;
//            for (k = 1; k < i-(j-bw+1)+1; k++)    //  correct, but w/ negative values
            for (k = 1; k < i-imin+1; k++)          //  should be correct w/o negative values
            {
//printf("\tG[%d][%d] += G[%d][%d]*G[%d][%d]\n", i,j, i-k,i, i-k,j);
//printf("\tG[%d][%d] += G[%d][%d]*G[%d][%d]\n", j-i,0, j-(i-k),0, j-(i-k),j-i);
                ods += G[j-(i-k)][0]*G[j-(i-k)][j-i];
            }
//printf("\tG[%d][%d] = (A[%d]-sum)/G[%d][%d]\n", i,j, j-i, i,i);
//printf("\tG[%d][%d] = (A[%d]-sum)/G[%d][%d]\n", j-i,0, j-i, j-i,j-i);
            G[j-i][0] = (A[j-i] - ods) / G[j-i][j-i];
            ds += G[j-i][0]*G[j-i][0];
        }
//printf("G[%d][%d] = A[%d]-sum\n", j,j, 0);
//printf("G[%d][%d] = A[%d]-sum\n", 0,0, 0);
assert (ds < A[0]);
        G[0][0] = sqrt(A[0] - ds);

        //  Convergence Test and save off pr
        norm = 0.0;
        for (i = 0; i < bw; i++)
        {
            norm += (G[i][0] - pr[i])*(G[i][0] - pr[i]);
            pr[i] = G[i][0];
        }
        norm = sqrt(norm);
        if (norm <= tol)
        {
printf("itr = %02d, norm = %e\n", j+1, norm);
            break;
        }

        //  Save off "previous row" (i.e. newly computed row)

//display_dynarr_d(G, bw, bw);
        //  Shift moving window
        for (i = bw-1; i > 0; i--)
        {
            for (k = i; k > 0; k--)
            {
                G[i][k] = G[i-1][k-1];
            }
        }
/*
        if (c[0] != 'r')
        {
            printf("Please press a key to step or \"r\" to run and then enter: ");
            scanf("%s", c);
        }
*/
    }

/*
    Solve (G^T)y = b using backward substitution
    Solve (G)x = y  using forward substitution
*/

    //  Free dynamically allocated memory
    dynfree(G[0]);
    dynfree(G);
    dynfree(pr);
}

// End of file