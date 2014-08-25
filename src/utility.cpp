//-------|---------|---------|---------|---------|---------|---------|---------|
/*
utility functions
*/

#include "stdafx.h"

/*
NOTE: https://software.intel.com/en-us/node/433552 for data alignment
NOTE: https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/
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

// End of file