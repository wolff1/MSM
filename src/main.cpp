//-------|---------|---------|---------|---------|---------|---------|---------|
/*
Defines the entry point for the console application.
*/

#include "stdafx.h"
#include "phiC1.h"	// Remove eventually

void test_mkl_MMM(void);
void mkl_memory_check(void);
void test_sinc(void);

int main(int argc, char* argv[])
{
	int choice = 1;

	while (choice > 0)
	{
		printf("************** MENU **************\n");
		printf("* 1 - Display Phi Test           *\n");
		printf("* 2 - Test MKL MMM               *\n");
		printf("* 3 - Display Gamma Tests        *\n");
		printf("* 4 - Produce Figure 1           *\n");
		printf("* 5 - Print nesting coefficients *\n");
		printf("* 6 - Produce Figure 4 (B-spline)*\n");
		printf("* 7 - Produce Figure 5 (C1)      *\n");
		printf("* 8 - Test theta and thetap      *\n");
		printf("* 9 - Test blurring operator     *\n");
		printf("*10 - Test polynomial multiply   *\n");
		printf("*11 - Test omega' values         *\n");
		printf("*12 - Test operator conversion   *\n");
		printf("*13 - Test omega values          *\n");
		printf("*14 - Produce Figure 2 (sinc)    *\n");
		printf("**********************************\n");
		printf("* 0 - Exit                       *\n");
		printf("**********************************\n");
		printf("Your selection: ");
		scanf("%d", &choice);

		// Give them what they want.
		switch (choice)
		{
			case 1:	// Test phi and phi'
				phi_test_all();
				break;

			case 2:	// Test matrix multiplication in MKL library
				test_mkl_MMM();
				break;

			case 3:	// Test each smoothing function independent of scaling
				gamma_test_all();
				break;

			case 4:	// This will produce plot (use OSX) for smoothing
				splitting_test();
				break;

			case 5:	// This will produce numerical values for J_n
				print_nesting_coefficients();
				break;

			case 6:	// This will produce plot (use OSX) for nesting
				phi_nesting_test();
				break;

			case 7:	// This will produce plots for C1 testing
				driverC1();
				break;

			case 8:	// test theta functions
				test_thetas();
				break;

			case 9:	// test blurring operator
				test_blurring_operator();
				break;

			case 10:	// test polynomial multiplication
				test_mpoly();
				break;

			case 11:	// test omega' values
				test_omegap();
				break;

			case 12:	// test operator basis conversion
				test_convert_to_shifts();
				break;

			case 13:	//	test omega values
				test_omega();
				break;

			case 14:
				test_sinc();
				break;

			case 0:	// Exit
				break;

			default:	// Invalid
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

/*
Compare cubic, quintic, and septic splines to sinc function
*/
void test_sinc(void)
{
	short		i = 0;
	short		j = 0;
	short		n = 0;
	short		p = 0;
	double		umin = 0.0;
	double		umax = 0.0;
	double		u = 0.0;
	short		samples = 100;
	short		omega_max = 100;
	double*		omega3 = NULL;
	double*		omega5 = NULL;
	double*		omega7 = NULL;
	double**	results = NULL;

	printf("Enter the number of samples: ");
	scanf("%hd", &samples);

	printf("Enter the number of omega values to compute: ");
	scanf("%hd", &omega_max);
	
	printf("Enter the minimum u value: ");
	scanf("%lf", &umin);

	printf("Enter the maximum u value: ");
	scanf("%lf", &umax);

	//	Dynamically allocate memory
	omega3 = (double*) dynvec(omega_max,sizeof(double));
	omega5 = (double*) dynvec(omega_max,sizeof(double));
	omega7 = (double*) dynvec(omega_max,sizeof(double));
	results = (double**) dynarr_d(4,samples+1);

	//	Compute omega values
	compute_omega(4, omega_max, omega3);
	compute_omega(6, omega_max, omega5);
	compute_omega(8, omega_max, omega7);

	for (i = 0; i <= samples; i++)
	{
		u = umin + (umax-umin)*(double)i /samples;

		//	Cubic
		p = 4;
		for (j = 0; j < p; j++)
		{
			n = (short)floor(u) - p/2 + 1 + j;
			assert(abs(n) < omega_max);
			results[0][i] += omega3[abs(n)]*phi(p,u-(double)n,NULL);
		}

		//	Quintic
		p = 6;
		for (j = 0; j < p; j++)
		{
			n = (short)floor(u) - p/2 + 1 + j;
			assert(abs(n) < omega_max);
			results[1][i] += omega5[abs(n)]*phi(p,u-(double)n,NULL);
		}

		//	Septic
		p = 8;
		for (j = 0; j < p; j++)
		{
			n = (short)floor(u) - p/2 + 1 + j;
			assert(abs(n) < omega_max);
			results[2][i] += omega7[abs(n)]*phi(p,u-(double)n,NULL);
		}

		//	Sinc function
		if (u != 0.0)
		{
			results[3][i] = sin(u*PI) / (u*PI);
		}
		else
		{
			results[3][i] = 1.0;
		}

		printf("u = %f, %f, %f, %f, %f\n", u, results[0][i], results[1][i], results[2][i], results[3][i]);
	}

	//	Free dynamically allocated memory
	dynfree(omega3);
	dynfree(omega5);
	dynfree(omega7);
	dynfree(results[0]);
	dynfree(results);
}

// End of file