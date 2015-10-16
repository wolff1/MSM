//-------|---------|---------|---------|---------|---------|---------|---------|
/*
memory.c - utility functions
*/

#include "memory.h"

/*
NOTE: https://software.intel.com/en-us/node/433552 for data alignment
NOTE: https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/
		GUID-637284D3-4D1F-4D6C-9509-382CB2DD1A3D.htm for routines involving
		linear systems
*/

/*
Dynamically allocate memory and zero it out
*/
void* dynmem(size_t size)
{
	//void*	ptr = calloc(size, sizeof(...));
	void*	ptr = mkl_malloc(size, MEM_ALIGN);
	if (size > 0)
	{	//	FIXME: Should assert get removed for performance once code is stable?
		assert(ptr != NULL);
	}
	memset(ptr, 0, size);
	return ptr;
}

/*
Free allocated memory
*/
void dynfree(void* ptr)
{
	//free(ptr);
	mkl_free(ptr);
}

/*
Use internal MKL functions to check for a memory leak
*/
void mkl_memory_check(void)
{
	MKL_INT64	AllocatedBytes;
	int			N_AllocatedBuffers;

	/**************************************************************************/
	// I'm not entirely sure which buffers this is freeing
	// This should only be called after last MKL usage
	mkl_free_buffers();

	AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
	if (AllocatedBytes > 0)
	{
		printf("MKL memory leak!\n");
		printf("After mkl_free_buffers there are %ld bytes in %d buffers\n",
			(long) AllocatedBytes, N_AllocatedBuffers);
	}

}

/*
Create 2D array of type double
*/
double** dynarr_d(unsigned long rows, unsigned long cols)
{
	unsigned long i = 0;
	double** ptr = NULL;
	double* buf = NULL;

	// Allocate pointer buffer
	ptr = (double**) dynvec(rows, sizeof(double));

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
Create vector (1D array) whose elements have size "size"
*/
void* dynvec(unsigned long rows, size_t size)
{
	return dynmem(rows*size);
}

/*
Display elements of 2D array
*/
void display_dynarr_d(double** A, unsigned long rows, unsigned long cols)
{
	unsigned long i = 0;
	unsigned long j = 0;

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
void display_vector_d(double* x, unsigned long rows)
{
	unsigned long i = 0;

	for (i = 0; i < rows; i++)
	{
		//printf("%f\n", x[i]);
		printf("%04lu\t%+e\n", i, x[i]);
	}
	printf("\n");
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
void bibst_lss(long max_itr, double tol,
				short A_len, double* A,
				short b_len, short b_nnz, double* b,
				short x_len, double* x)
{
	lapack_int			i = 0;
	lapack_int			j = 0;
	char				uplo = 'L';
	lapack_int			m = 4;					//	Experimentally determined multiplier to find minimum required length to still attain machine precision in solution
	lapack_int			n = MAX(m*b_len,x_len);	//	Order of A [~ 4*len(b_len)]
	lapack_int			kd = A_len-1;			//	# of super OR sub diagonals in A [= A_len - 1]
	double**			ab = NULL;				//	(kd+1)*n matrix represented banded A
	lapack_int			ldab = kd+1;			//	Leading dimension of ab (ldab >= kd + 1)
	lapack_int			info = 0;				//	Return code from lapack routine
	lapack_int			nrhs = 1;				//	# of right-hand sides
	double*				rhs = NULL;				//	Right-hand side vector used in solve

	//	Allocate memory for banded storage of A
	ab = (double**) dynarr_d(ldab,2*n+1);

	//	Create banded representation of A for MKL (only lower portion of A, from diagonal to furthest sub-diagonal)
	for (i = 0; i < ldab; i++)
	{
		for (j = 0; j < 2*n+1; j++)
		{
			ab[i][j] = A[i];
		}
	}
//display_dynarr_d(ab,ldab,n);

	//	Compute Cholesky factors of symmetric, real, banded, positive-definite matrix A
	if (!(info = LAPACKE_dpbtrf(LAPACK_ROW_MAJOR, uplo, 2*n+1, kd, ab[0], 2*n+1)))
	{
//display_dynarr_d(ab,ldab,n);

		//	Allocate memory for rhs vector
		rhs = (double*) dynvec(2*n+1, sizeof(double));

		//	Fill in rhs vector with values from b vector
//display_vector_d(b, b_len);

		rhs[n] = b[0];
		for (i = 1; i < MIN(b_nnz,n); i++)
		{
			rhs[n+i] = b[i];
			rhs[n-i] = b[i];
		}
//display_vector_d(rhs, 2*n+1);

		//	Use A and b to solve for x
		if (!(info = LAPACKE_dpbtrs(LAPACK_ROW_MAJOR, uplo, 2*n+1, kd, nrhs, ab[0], 2*n+1, rhs, (lapack_int)1)))
		{
//display_vector_d(rhs, 2*n+1);

			//	Success -- fill b into x (or something)
			for (i = 0; i < x_len; i++)
			{
				x[i] = rhs[n+i];
			}
//display_vector_d(x, x_len);
		}

		//	Free dynamically allocated memory for rhs vector
		dynfree(rhs);
	}

	//	Free dynamically allocated memory
	dynfree(ab[0]);
	dynfree(ab);
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
void bibst_lss2(long max_itr, double tol,
				short A_len, double* A,
				short b_len, short b_nnz, double* b,
				short x_len, double* x)
{
    short           i = 0;
    short           j = 0;
    short           k = 0;
    short           bw = A_len; //  Change A_len to bw in definition?
    double          ods = 0.0;  //  off-diagonal sum
    double          ds = 0.0;   //  on-diagonal sum
    double          norm = 0.0; //  Used in convergence check
    double**        G = NULL;   //  Lower triangular Cholesky factor
    double*         pr = NULL;  //  Previous Row (used in infinite routine)
	double*			y = NULL;	//	Used in backward substition
    double          A1 = 0.0;
    
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
    y = (double*) dynvec(MAX(b_nnz,x_len), sizeof(double));

    if (max_itr < 0)
        max_itr = 1000;

/*
    Factor A = (G^T)(G) (Upper-Lower Cholesky)
        -> Use infinite Cholesky algorithm to get "converged" row
        -> Use finite Cholesky algorithm to get "unique" rows
*/
    //  "INFINITE" CHOLESKY
    for (j = 0; j < max_itr; j++)
    {
        ds = 0.0;	//ds -> diagonal sum: sum used in computing the diagonal
        for (i = MIN(j,bw-1); i > 0; i--)//  For j < bw-1, we want A[i] = 0.0
        {
            ods = 0.0;	//ods -> off-diagonal sum, used in off-diagonal values
            for (k = 1; k < bw-i; k++)
            {
                ods += G[i+k][0]*G[i+k][i];
            }
            G[i][0] = (A[i] - ods) / G[i][i];
            ds += G[i][0]*G[i][0];
        }
        G[0][0] = sqrt(A[0] - ds);

        //  Convergence Test & save off "previous row" (newly computed row)
        norm = 0.0;
        for (i = 0; i < bw; i++)
        {
            norm += (G[i][0] - pr[i])*(G[i][0] - pr[i]);
            pr[i] = G[i][0];
        }
        norm = sqrt(norm);
        if (norm <= tol)
        {
            break;
        }

        //  Shift moving window
        for (i = bw-1; i > 0; i--)
        {
            for (k = i; k > 0; k--)
            {
                G[i][k] = G[i-1][k-1];
            }
        }
    }

    //  Reset "infinite" Cholesky factor to be all converged row
    for (i = 0; i < bw; i++)
    {
        for (j = 0; j <= i; j++)
        {
            G[i][j] = pr[i-j];
        }
    }

    //  FINITE CHOLESKY
    for (j = ceil(bw/2.0)-1; j >= 0; j--)    //  Loop over remaining columns
    {
        ds = 0.0;	// diagonal sum - used in computed diagonal value
        for (i = bw-1; i > 0; i--)
        {
            ods = 0.0;	// off-diagonal sum, used for off-diagonal values
            for (k = 1; k < bw-i; k++)
            {
				ods += G[i+k][0]*G[i+k][i];
            }

            //  Determine true value derived from A to use
            if ((j > 0) && (i < 2*(ceil(bw/2.0)-j-bw%2) + bw%2))
            {
                A1 = A[i] + A[(bw-1)-abs(bw-2*j-i-1)];
            }
            else
            {
                A1 = A[i];
            }
            G[i][0] = (A1 - ods) / G[i][i];
            ds += G[i][0]*G[i][0];
        }

        //  Determine true value derived from A to use
        if (j > 0)
        {
            //  NOTE: Same formula as above, but i = 0
            A1 = A[0] + A[(bw-1)-abs(bw-2*j-1)];
        }
        else
        {
            A1 = A[0]/2.0;
        }
        G[0][0] = sqrt(A1 - ds);

		if (j > 0)
		{
			//  Shift moving window
			for (i = bw-1; i > 0; i--)
			{
				for (k = i; k > 0; k--)
				{
					G[i][k] = G[i-1][k-1];
				}
			}
		}
    }

/*
    Solve (G^T)y = b using backward substitution
    Solve (G)x = y  using forward substitution
*/
	if (b_nnz > b_len)
	{	//	Sanity check - avoid stepping out of bounds
		b_nnz = b_len;
	}

	//	BACKWARD SUBSTIUTION
	for (i = b_nnz-1; i >= bw; i--)
	{
		//  NOTE:   First divide b by 2 b/c A was implicitely divided by 2
		y[i] = 0.5*b[i];
		for (j = 1; j < MIN(b_nnz-i,bw); j++)
		{	// pr is converged row in Cholesky factor
			y[i] -= pr[j]*y[i+j];
		}
		y[i] /= pr[0];	//	pr[0] is diagonal of the converged row
	}

	k = MIN(bw,b_nnz);
	for (i = k-1; i >= 0; i--)
	{
		y[i] = 0.5*b[i];
		for (j = 1; j < k-i; j++)
		{
			y[i] -= G[j][i]*y[i+j];
		}

		// NOTE: MIN() keeps j from exceeding b_nnz-th column
		for (j = k-i; j < MIN(b_nnz-i,bw); j++)
		{
			y[i] -= pr[j]*y[i+j];
		}

		y[i] /= G[i][i];
	}

	//	FORWARD SUBSTITUION
	k = MIN(bw,x_len);
	for (i = 0; i < k; i++)
	{
		x[i] = y[i];
		for (j = i; j > 0; j--)
		{
		  x[i] -= G[i][i-j]*x[i-j];
		}
		x[i] /= G[i][i];
	}

	for (i = k; i < x_len; i++)
	{
		x[i] = y[i];
		for (j = bw-1; j > 0; j--)
		{
		  x[i] -= pr[j]*x[i-j];
		}
		x[i] /= pr[0];
	}

	//  Free dynamically allocated memory
    dynfree(G[0]);
    dynfree(G);
    dynfree(pr);
	dynfree(y);
}

// End of file
