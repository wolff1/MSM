//-------|---------|---------|---------|---------|---------|---------|---------|
/*
Defines the entry point for the console application.
*/

#include "stdafx.h"

void splitting_test(void);
void plot_splittings(int samples, int nlev, int k, double a, double d, double* X, double** F);

int main(int argc, char* argv[])
{
/*
	int			i = 0;
	int			samples = 0;
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	char*		data_file = NULL;
	int			p = 0;		// Phi
	int			k = 0;		// Gamma
	double*		c = NULL;	// Gamma
	//// MKL
	//double**	A = NULL;
	//double**	B = NULL;
	//double**	C = NULL;
*/
	MKL_INT64	AllocatedBytes;
	int			N_AllocatedBuffers;

	splitting_test();

#if 0
	// Test plotting routine and gnuplot
//	plotf("damped_sine_commands.txt", "damped_sine_data.txt");

	/****************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%d", &p);
	assert(p % 2 == 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%d", &k);
	assert(k > 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));
	c = (double*) dynvec(k+1,sizeof(double));
	data_file = (char*) dynvec(strlen(GP_DATA_DIR) +
								MAX(strlen(PHI_DATA),strlen(GAMMA_DATA)) + 1,
								sizeof(char));
	/****************************************************************************/
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
	sprintf(data_file, "%s%s", GP_DATA_DIR, PHI_DATA);
	printf("Plotting Phi(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plots2d(samples, X, F, data_file);
	//	Show Phi'(x)
	printf("Plotting Phi'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	/****************************************************************************/
	// Test softener
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
	//	Show Gamma(x)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GAMMA_DATA);
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plots2d(samples, X, F, data_file);
	//	Show Gamma'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	/****************************************************************************/
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

	/****************************************************************************/
	// Test Theta* - short range splitting function
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = theta(c, k, X[i], &DF[i], 0);
/*
		printf("%02d:\ttheta*(%f) = %f\ttheta*'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Theta*(x)
	sprintf(data_file, "%s%s%s", GP_DATA_DIR, "star_", THETA_DATA);
	printf("Plotting theta*(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plots2d(samples-1, &X[1], &F[1], data_file);
	//	Show Theta*'(x)
	printf("Plotting theta*'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	// Test Theta - long range splitting function (finite)
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = theta(c, k, X[i], &DF[i], 1);
/*
		printf("%02d:\ttheta(%f) = %f\ttheta'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Theta(x)
	sprintf(data_file, "%s%s", GP_DATA_DIR, THETA_DATA);
	printf("Plotting theta(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plots2d(samples-1, &X[1], &F[1], data_file);
	//	Show Theta'(x)
	printf("Plotting theta'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	// Test Theta - long range splitting function (infinite)
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = gamma(c, k, X[i], &DF[i]);
/*
		printf("%02d:\tgamma(%f) = %f\tgamma'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show gamma(x)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GAMMA_DATA);
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plots2d(samples-1, &X[1], &F[1], data_file);
	//	Show gamma'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	/****************************************************************************/
	// Free allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F);
	dynfree(DF);
	dynfree(data_file);
#endif

	/****************************************************************************/
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

void splitting_test(void)
{
	int			i = 0;
	int			samples = 0;
	double*		c = NULL;	// Gamma
	int			k = 0;
	double		a = 0.0;
	double		one_over_a = 0.0;
	double		a1 = 0.0;
	int			nlev = 0;
	int			l = 0;
	double		d = 0.0;

	double*		X = NULL;
	double**	F = NULL;
	double**	DF = NULL;
	char*		data_file = NULL;

	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter leves parameter nlev = ");
	scanf("%d", &nlev);
	assert(nlev > 1);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%d", &k);
	assert(k > 0);

	// Get a where a is cut-off distance
	printf("Please enter cut-off parameter a = ");
	scanf("%f", &a);
	assert(a > 0.0);
	one_over_a = 1.0/a;

	// Get d where d is domain
	printf("Please enter domain parameter d = ");
	scanf("%f", &d);
	assert(d > 0.0);

	// Dynamically allocate memory
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double**) dynarr_d(nlev+1,samples+1);
	DF = (double**) dynarr_d(nlev+1,samples+1);
	c = (double*) dynvec(k+1,sizeof(double));
	data_file = (char*) dynvec(strlen(GP_DATA_DIR) + MAX(strlen(PHI_DATA),strlen(GAMMA_DATA)) + 1, sizeof(char));

	// Initialize gamma coefficients
	gamma_init(k, c);

	// Evaluate different splitting functions over domain
	for (i = 1; i <= samples; i++)
	{
		X[i] = (d*(double)i/(double)samples);
		// Theta* - Short range part of splitting (finite)
		a1 = one_over_a;
		F[0][i] = theta(c, k, a1*X[i], &DF[0][i], 0);
		F[0][i] = a1*F[0][i];
		DF[0][i] = a1*a1*DF[0][i];
		for (l = 1; l < nlev - 1; l++)
		{
			// Theta - Intermediate long range part(s) of splitting (finite)
			a1 = a1*one_over_a;
			F[l][i] = theta(c, k, a1*X[i], &DF[l][i], 1);
			F[l][i] = a1*F[l][i];
			DF[l][i] = a1*a1*DF[l][i];
		}
		// Gamma - Top level long range part of splitting (infinit)
		a1 = a1*one_over_a;
		F[l][i] = gamma(c, k, a1*X[i], &DF[l][i]);
		F[l][i] = a1*F[l][i];
		DF[l][i] = a1*a1*DF[l][i];
		// Kernel: 1/X and Kernel' = -1/X^2
		F[nlev][i] = 1.0/X[i];
		DF[nlev][i] = -F[nlev][i]*F[nlev][i];
	}

	// Plot splittings on single graph along with f(x) = 1/x
	plot_splittings(samples, nlev, k, a, d, X, F);
	plot_splittings(samples, nlev, k, a, d, X, DF);

	// Free dynamically allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F[0]);
	dynfree(F);
	dynfree(DF[0]);
	dynfree(DF);
	dynfree(data_file);
}

/*
samples:	number of data points
nlev:		number of levels (3->2 grids)
k:			continuity of smoothing
a:			cut-off distance
d:			length of domain starting at 0.0
X:			independent variable, 0.0 <= X <= d (vector, 1 x samples)
F:			dependent variable(s) (array, nlev x samples)
*/
void plot_splittings(int samples, int nlev, int k, double a, double d, double* X, double** F)
{
	FILE*	data = NULL;
	FILE*	cmd = NULL;

	// Create DATA file with X, F[0], F[1], ..., F[nlev] columns and samples rows
	data = fopen("data_file_name", "w");
	assert(data != NULL);
//FIXME
	// Create COMMAND file to plot all of the columns above
	cmd = fopen("cmd_file_name", "w");
	assert(cmd != NULL);
//FIXME

	// Call gnuplot to plot DATA file using COMMAND file
	plotf2d("cmd_file_name", "data_file_name");

	if (fclose(data))
	{
		printf("Error closing DATA file.\n");
	}

	if (fclose(cmd))
	{
		printf("Error closing DATA file.\n");
	}

	// Delete data file
	if (remove("data_file_name"))
	{
		printf("Error removing DATA file <%s>.\n", "data_file_name");
	}

	// Delete command data file
	if (remove("cmd_file_name"))
	{
		printf("Error removing COMMAND file <%s>.\n", "cmd_file_name");
	}
}

/*
#define GP_DATA_DIR			"../../../data/"
#define GP_CMD_DIR			"../../../gnuplot/"
#define GP_DATA_TMP			"tmp_data.dat"
#define GP_CMD_TEMPLATE		"tmp_data_template_win.txt"
*/
// End of file