//-------|---------|---------|---------|---------|---------|---------|---------|
/*
Defines the entry point for the console application.
*/

#include "stdafx.h"

void test_mkl_MMM(void);
void mkl_memory_check(void);

int main(int argc, char* argv[])
{
	int choice = 1;

	while (choice > 0)
	{
		printf("********** MENU ***********\n");
		printf("* 1 - Display Phi Test    *\n");
		printf("* 2 - Test MKL MMM        *\n");
		printf("* 3 - Display Gamma Tests *\n");
		printf("* 4 - Produce Figure 1    *\n");
		printf("* 0 - Exit                *\n");
		printf("***************************\n");
		printf("Your selection: ");
		scanf("%d", &choice);

		// Give them what they want.
		switch (choice)
		{
			case 1:
				// Test phi and phi'
				phi_test_all();
				break;

			case 2:
				// Test matrix multiplication in MKL library
				test_mkl_MMM();
				break;

			case 3:
				// Test each smoothing function independent of method scaling
				gamma_test_all();
				break;

			case 4:
				// This will produce plot (use OSX) for smoothing
				splitting_test();
				break;

			case 0:
				break;

			default:
				printf("You selected <%d> which is INVALID\n", choice);
				choice = 1;
		}

		if (choice > 0)
		{
			// Put some space before next menu is displayed
			printf("\n\n");
		}
	}

	// Use MKL functions to check for memory leak
	mkl_memory_check();

	return 0;
}

/*
Use MKL and CBLAS to do matrix-matrix multiplication
*/
void test_mkl_MMM(void)
{
	MKL_INT64	AllocatedBytes;
	int			N_AllocatedBuffers;
	double**	A = NULL;
	double**	B = NULL;
	double**	C = NULL;

	/**************************************************************************/
	// Test Matrix-Matrix multiplication (using MKL)
	printf("\n");
	A = dynarr_d(2,2);
	A[0][0] = 2.0;	A[0][1] = 0.0;
	A[1][0] = 0.0;	A[1][1] = 4.0;
	B = dynarr_d(2,2);
	B[0][0] = 0.5;	B[0][1] = 0.0;
	B[1][0] = 0.0;	B[1][1] = 0.25;
	C = dynarr_d(2,2);
	C[0][0] = 0.0;	C[0][1] = 0.0;
	C[1][0] = 0.0;	C[1][1] = 0.0;

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0,
				A[0], 2, B[0], 2, 0.0, C[0], 2);

	display_dynarr_d(A,2,2);
	display_dynarr_d(B,2,2);
	display_dynarr_d(C,2,2);

	AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
	printf("DGEMM uses %ld bytes in %d buffers\n",
			(long) AllocatedBytes, N_AllocatedBuffers);

	dynfree(A[0]);
	dynfree(A);
	dynfree(B[0]);
	dynfree(B);
	dynfree(C[0]);
	dynfree(C);
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

// End of file