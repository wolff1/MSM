//-------|---------|---------|---------|---------|---------|---------|---------|
/*
Defines the entry point for the console application.
*/

#include "stdafx.h"

void splitting_test(void);
void plot_splittings(int samples, int nlev, int k, double a, double d, double* X, double** F);

int main(int argc, char* argv[])
{
	MKL_INT64	AllocatedBytes;
	int			N_AllocatedBuffers;

	splitting_test();
#if 0
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



	// Test plotting routine and gnuplot
//	plotf("damped_sine_commands.txt", "damped_sine_data.txt");

	/**************************************************************************/
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
	data_file = (char*) dynvec(GP_DATA_DIR_LEN +
							   MAX(PHI_DATA_len,GAMMA_DATA_len) + 1,
							   sizeof(char));
	/**************************************************************************/
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

	/**************************************************************************/
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

	/**************************************************************************/
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

	/**************************************************************************/
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

	/**************************************************************************/
	// Free allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F);
	dynfree(DF);
	dynfree(data_file);
#endif

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

	double		f = 0.0;
	double		df = 0.0;
	double		tol = pow(0.1,13);

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
	scanf("%lf", &a);
	assert(a > 0.0);
	one_over_a = 1.0/a;

	// Get d where d is domain
	printf("Please enter domain parameter d = ");
	scanf("%lf", &d);
	assert(d > 0.0);

	// Dynamically allocate memory
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double**) dynarr_d(nlev+1,samples+1);
	DF = (double**) dynarr_d(nlev+1,samples+1);
	c = (double*) dynvec(k+1,sizeof(double));
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + MAX(PHI_DATA_LEN,GAMMA_DATA_LEN) + 1, sizeof(char));

	// Initialize gamma coefficients
	gamma_init(k, c);

	// Evaluate different splitting functions over domain
	for (i = 0; i <= samples; i++)
	{
f = 0.0;	// sanity check
df = 0.0;	// sanity check
		X[i] = (d*(double)i/(double)samples);
		// Theta* - Short range part of splitting (finite)
		a1 = one_over_a;
		F[0][i] = theta(c, k, a1*X[i], &DF[0][i], 0);
		F[0][i] = a1*F[0][i];
		DF[0][i] = a1*a1*DF[0][i];
f += F[0][i];
df += DF[0][i];
		for (l = 1; l < nlev - 1; l++)
		{
			// Theta - Intermediate long range part(s) of splitting (finite)
			F[l][i] = theta(c, k, a1*X[i], &DF[l][i], 1);
			F[l][i] = a1*F[l][i];
			DF[l][i] = a1*a1*DF[l][i];
f += F[l][i];
df += DF[l][i];
			a1 = a1*one_over_a;
		}
		// Gamma - Top level long range part of splitting (infinit)
		a1 = a1*one_over_a;
		F[l][i] = gamma(c, k, a1*X[i], &DF[l][i]);
		F[l][i] = a1*F[l][i];
		DF[l][i] = a1*a1*DF[l][i];
f += F[l][i];
df += DF[l][i];
		// Kernel: 1/X and Kernel' = -1/X^2
		F[nlev][i] = 1.0/X[i];
		DF[nlev][i] = -F[nlev][i]*F[nlev][i];
if (i > 0)
{
	if (fabs(F[nlev][i] - f) >= tol)
		printf("F = %f, f = %f, |.| = %e\n", F[nlev][i], f, fabs(F[nlev][i] - f));
	if (fabs(DF[nlev][i] - df) >= tol)
		printf("DF = %f, df = %f, |.| = %e\n", DF[nlev][i], df, fabs(DF[nlev][i] - df));
	assert(fabs(F[nlev][i] - f) < tol);
	assert(fabs(DF[nlev][i] - df) < tol);
}
	}

	// Plot splittings on single graph along with f(x) = 1/x
	plot_splittings(samples, nlev, k, a, d, X, F);
//	plot_splittings(samples, nlev, k, a, d, X, DF);

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
NOTE:		Ignore first sample becuase f(0) = 1/0 is undefined!
*/
void plot_splittings(int samples, int nlev, int k, double a, double d, double* X, double** F)
{
	FILE*	data = NULL;
	FILE*	cmd = NULL;
	char*	data_file = NULL;
	char*	cmd_file = NULL;
	char*	cmd_file_name = "splitting_";
	char*	cmd_file_extension = ".txt";
	int		cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN + strlen(cmd_file_extension);

	// File stuff
	char* buf = NULL;
	size_t bufmax = 0;
	char* buf2 = NULL;
	size_t buf2max = 0;
	size_t buflen = 0;
	size_t bytes = 0;
	int i = 0;
	int l = 0;

	assert(X != NULL);
	assert(F != NULL);

	// Dynamically allocate memory
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + GP_DATA_TMP_LEN + 1,sizeof(char));
	cmd_file = (char*) dynvec(GP_CMD_DIR_LEN + cmd_file_name_len + 1,sizeof(char));

	// Build file name(s)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(cmd_file, "%s%s%s%s", GP_CMD_DIR, cmd_file_name, GP_TERM, cmd_file_extension);

	// Create DATA file with X, F[0], F[1], ..., F[nlev] columns and samples rows
	data = fopen(data_file, "w");
	assert(data != NULL);

	// Write DATA file
	bufmax = 64*(nlev+2);	// nlev + 2 columns, each a max of 64 chars wide 
	buf = (char*) dynvec(bufmax+1,sizeof(char));	// + 1 for NULL
	// Adjust F[0][0] and F[nlev][0] to be large number
	F[0][0] = 1000.0;
	F[nlev][0] = 1000.0;
	for (i = 0; i <= samples; i++)
	{
// FIXME - There is probably a better way to format this file... binary data would be most accurate, right?
		// Add X to buffer
		buflen = sprintf(buf, "%.32f", X[i]);
		assert(bufmax - buflen > 0);

		for (l = 0; l < nlev; l++)
		{
			// Add F[l][i] to buffer
			buflen += sprintf(&buf[buflen], " %.32f", F[l][i]);
			assert(bufmax - buflen > 0);
		}
		// Add F[nlev][i] to buffer
		buflen += sprintf(&buf[buflen], " %.32f\n", F[nlev][i]);
		assert(bufmax - buflen > 0);

		// Write buffer to file
		bytes = fwrite(buf, sizeof(char), buflen, data);
		if (bytes < buflen)
		{
			printf("<%d> bytes written to temporary file <%s> instead of <%d>\n", bytes, data_file, buflen);
			break;
		}

		// Clear out buffer for next usage
		memset(buf, 0, bufmax+1);
	}

	// Close DATA file
	if (fclose(data))
	{
		printf("Error closing DATA file.\n");
	}

	// Create COMMAND file to plot all of the columns above
	cmd = fopen(cmd_file, "w");
	assert(cmd != NULL);

	// Write CMD file
	buf2max = 256 + 32*nlev;
	buf2 = (char*) dynvec(buf2max+1,sizeof(char));
	// Create COMMAND file buffer
	buflen = sprintf(buf2,
		"set term %s\n"
		"set xlabel 'x'\n"
		"set ylabel 'f(x)'\n"
		"set title 'Splitting for %d-level MSM'\n"
		"set grid\n"
		"set style data lines\n"
		"set yrange [ -1.0 : %d ]\n"
		"plot "
		, GP_TERM, nlev, 10);//ceil(1.5*F[1][1]));
printf("|buf2| = %d\n", buflen);
	assert(buf2max - buflen > 0);

	// Plot a line for each level (columns 1 -> nlev+1)
	for (i = 0; i < nlev; i++)
	{
		buflen += sprintf(&buf2[buflen],
			"data_file using 1:%d,",
			i+2);
printf("|buf2| = %d\n", buflen);
		assert(buf2max - buflen > 0);
	}

	// Plot line for 1/r (column nlev+2) and finish file
	buflen += sprintf(&buf2[buflen],
		"data_file using 1:%d\n"
		"pause -1\n"
		"quit\n",
		nlev+2);
printf("|buf2| = %d\n", buflen);
	assert(buf2max - buflen > 0);

	// Write buffer to file
	bytes = fwrite(buf2, sizeof(char), buflen, cmd);
	if (bytes < buflen)
	{
		printf("<%d> bytes written to temporary file <%s> instead of <%d>\n", bytes, cmd_file, buflen);
	}

	// Close CMD file
	if (fclose(cmd))
	{
		printf("Error closing DATA file.\n");
	}

	// Call gnuplot to plot DATA file using COMMAND file
	plotf2d(cmd_file, data_file);

	// Delete data file
	if (remove(data_file))
	{
		printf("Error removing DATA file <%s>.\n", data_file);
	}

/*
	// Delete command data file
	if (remove(cmd_file))
	{
		printf("Error removing COMMAND file <%s>.\n", cmd_file);
	}
*/
	// Free dynamically allocated memory
	dynfree(data_file);
	dynfree(cmd_file);
	dynfree(buf);
	dynfree(buf2);
}

// End of file