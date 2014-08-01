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

// End of file