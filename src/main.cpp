//-------|---------|---------|---------|---------|---------|---------|---------|
/*
Defines the entry point for the console application.
*/

#include "stdafx.h"

//int _tmain(int argc, _TCHAR* argv[]) // WINDOWS SPECIFIC
int main(int argc, char* argv[])
{
	int i = 0;
	int samples = 0;
	double* X = NULL;
	double* F = NULL;
	double* DF = NULL;
	int p = 0;			// Phi
	int k = 0;			// Gamma
	double* c = NULL;	// Gamma
	//// MKL
	//double** A = NULL;
	//double** B = NULL;
	//double** C = NULL;
	MKL_INT64 AllocatedBytes;
	int N_AllocatedBuffers;

	// Test plotting routine and gnuplot
//	plotf("damped_sine_commands.txt", "damped_sine_data.txt");

	// Get number of data points to record
	printf("Please enter the number of samples = ");
//	scanf_s("%d", &samples, 1);
	scanf("%d", &samples);
	assert(samples > 0);

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
//	scanf_s("%d", &p, 1);
	scanf("%d", &p);
	assert(p % 2 == 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
//	scanf_s("%d", &k, 1);
	scanf("%d", &k);
	assert(k > 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));

	// Test centered B-spline where p must be even
	for (i = 0; i <= samples; i++)
	{
		X[i] = (double)(-p/2.0) + ((double)i/(double)samples)*((double)p);
		F[i] = phi(p, X[i], &DF[i]);
/*
		printf("%02d:\tphi(%f) = %f\tphi'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Phi(x)
	printf("Plotting Phi(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plots2d(samples, X, F, "phi.dat");
	//	Show Phi'(x)
	printf("Plotting Phi'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	// Test softener
	c = (double*) dynvec(k+1,sizeof(double));
	gamma_init(k, c);
	for (i = 0; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = gamma(c, k, X[i], &DF[i]);
/*
		printf("%02d:\tgamma(%f) = %f\tgamma'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Phi(x)
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plots2d(samples, X, F, "gamma.dat");
	//	Show Phi'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	//// Test Matrix-Matrix multiplication (using MKL)
	//printf("\n");
	//A = dynarr_d(2,2);
	//A[0][0] = 2.0;	A[0][1] = 0.0;
	//A[1][0] = 0.0;	A[1][1] = 4.0;
	//B = dynarr_d(2,2);
	//B[0][0] = 0.5;	B[0][1] = 0.0;
	//B[1][0] = 0.0;	B[1][1] = 0.25;
	//C = dynarr_d(2,2);
	//C[0][0] = 0.0;	C[0][1] = 0.0;
	//C[1][0] = 0.0;	C[1][1] = 0.0;

	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0,
	//			A[0], 2, B[0], 2, 0.0, C[0], 2);

	//display_dynarr_d(A,2,2);
	//display_dynarr_d(B,2,2);
	//display_dynarr_d(C,2,2);

	//AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
	//printf("DGEMM uses %ld bytes in %d buffers\n",
	//		(long) AllocatedBytes, N_AllocatedBuffers);

	//dynfree(A[0]);
	//dynfree(A);
	//dynfree(B[0]);
	//dynfree(B);
	//dynfree(C[0]);
	//dynfree(C);

	// Free allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F);
	dynfree(DF);

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

	return 0;
}

// End of file