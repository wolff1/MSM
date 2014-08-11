//-------|---------|---------|---------|---------|---------|---------|---------|
/*
gamma.cpp - the softening function used to split the kernel
*/

#include "stdafx.h"

/*
Solves for coefficients of gamma
k is degree of continuity
b will return k+1 vector of coefficients
*/
void gamma_init(int k, double* x)
{
	double** A = NULL;	// Coefficient matrix
	double* b = x;		// Alias output parameter x
	int i = 0;
	int j = 0;
	lapack_int rc = 0;
	lapack_int* piv = NULL;

	// b will first act as RHS, then return coefficients to caller
	assert(b != NULL);

	// Allocate memory for A
	A = dynarr_d(k+1, k+1);

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
/*
	// Display A
	printf("A = \n");
	display_dynarr_d(A,k+1, k+1);

	// Display b
	printf("b = \n");
	display_vector_d(b, k+1);
*/
	// Solve system Ax = b
	// NOTE: b is overwritten with x, A is overwritten with LU
	piv = (lapack_int*) dynvec(k+1,sizeof(lapack_int));
	rc = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) k+1, (lapack_int) 1,
						A[0], (lapack_int) k+1, piv, b, (lapack_int) 1);
	dynfree(piv);
	assert(rc == 0); // zero is SUCCESS
/*
	// Display c
	printf("c = \n");
	display_vector_d(b, k+1);
*/
	// Free allocated memory
	dynfree(A[0]);
	dynfree(A);
}

/*
Evaluate gamma and gamma' at (positive) position x
c is coefficient vector for gamma
*/
double gamma(double *c, int k, double x, double* dgamma)
{
	double f = 0.0;
	double df = 0.0;
	double xx = x*x;
	int i = 0;

	// gamma is a spline which becomes f(x) = 1/x for x >= 1.0
	if (x >= 1.0)
	{
		f = 1.0/x;
		df = -f*f;	//	-f*f = -1.0/(x*x);
	}
	else
	{
		// Use Horner's rule to evaluate polynomial
		f = c[k];			// Even powers
		df = c[k]*(2*k);	// Odd powers
		for (i = k-1; i >= 1; i--)
		{
			f = f*xx + c[i];
			df = df*xx + c[i]*(2*i);
		}
		f = f*xx + c[0];
		df = df*x;
	}

	if (dgamma != NULL)
		*dgamma = df;
	return f;
}

/*
Evaluate theta* and theta*' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta_star(double *c, int k, double x, double* dtheta_star)
{
    double f = 0.0;
    double df = 0.0;

    f = 1.0/x - gamma(c, k, x, &df);
    df = -1.0/(x*x) - df;

    if (dtheta_star != NULL)
        *dtheta_star = df;
    return f;
}

/*
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta(double *c, int k, double x, double* dtheta)
{
	double f = 0.0;
	double df = 0.0;
	double df2 = 0.0;

    f = gamma(c, k, x, &df) - 0.5*gamma(c, k, 0.5*x, &df2);
    df = df - 0.25*df2;

	if (dtheta != NULL)
		*dtheta = df;
	return f;
}

/*** DRIVER FUNCTIONS BELOW ***/

/*
Prompts user for samples, number of levels, degree of continuity for smoothing,
	cut-off value and domain length.
Then uses theta_star, theta, and gamma to compute (and plot) the smoothings for
	the method.
*/
void splitting_test(void)
{
	// Static memory variables
	int			i = 0;
	int			samples = 0;
	int			k = 0;
	double		a = 0.0;
	double		one_over_a = 0.0;
	double		al = 0.0;
	int			nlev = 0;
	int			l = 0;
	double		d = 0.0;
	double		f = 0.0;	//FIXME - Remove f, df, tol, and tol check below
	double		df = 0.0;
	double		tol = pow(0.1,15);
	// Dynamic memory variables
	double*		c = NULL;	// Gamma
	double*		X = NULL;
	double**	F = NULL;
	double**	DF = NULL;
	char*		data_file = NULL;

	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter levels parameter nlev = ");
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
		al = one_over_a;
		F[0][i] = al*theta_star(c, k, al*X[i], &DF[0][i]);
		DF[0][i] = al*al*DF[0][i];
f += F[0][i];
df += DF[0][i];

		for (l = 1; l < nlev - 1; l++)
		{
			// Theta - Intermediate long range part(s) of splitting (finite)
			F[l][i] = al*theta(c, k, al*X[i], &DF[l][i]);
			DF[l][i] = al*al*DF[l][i];
f += F[l][i];
df += DF[l][i];

			al = 0.5*al;
		}
		// Gamma - Top level long range part of splitting (infinit)
		F[l][i] = gamma(c, k, al*X[i], &DF[l][i]);
		F[l][i] = al*F[l][i];
		DF[l][i] = al*al*DF[l][i];
f += F[l][i];
df += DF[l][i];

		// Kernel: 1/X and Kernel' = -1/X^2
		F[nlev][i] = 1.0/X[i];
		DF[nlev][i] = -F[nlev][i]*F[nlev][i];
if (i > 0)
{
    double Frel = fabs(F[nlev][i] - f)/fabs(F[nlev][i]);
    double DFrel = fabs(DF[nlev][i] - df)/fabs(DF[nlev][i]);
	if (Frel >= tol)
		printf("i = %d, X = %f, F = %f, f = %f, |.| = %e\n", i, X[i], F[nlev][i], f, Frel);
	if (DFrel >= tol)
		printf("i = %d, X = %f, DF = %f, df = %f, |.| = %e\n", i, X[i], DF[nlev][i], df, DFrel);
	assert(Frel < tol);
	assert(DFrel < tol);
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
*/
void plot_splittings(int samples, int nlev, int k, double a, double d,
						double* X, double** F)
{
	// Static memory variables
	char*	cmd_file_name = "splitting_";
	char*	cmd_file_extension = ".txt";
	size_t	cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN + strlen(cmd_file_extension);
	// File stuff
	size_t  bufmax = 0;
	size_t  buf2max = 0;
	size_t  buflen = 0;
	size_t  bytes = 0;
	int		i = 0;
	int		l = 0;
	// Dynamic memory variables
	FILE*	data = NULL;
	FILE*	cmd = NULL;
	char*	data_file = NULL;
	char*	cmd_file = NULL;
	char*	buf = NULL;
	char*	buf2 = NULL;
    char*   buf2_1 =    "set term %s\n"
                        "set xlabel 'r'\n"
//                      "set ylabel 'g_l(r)'\n" // below two lines rotates by -90 degrees
                        "set lmargin 10\n"
                        "set label 1 'g_l(r)' at graph -0.1, graph 0.5\n"
                        "set label 2 '(a)' at graph 0.39, graph -0.08\n"
                        "set label 3 '(2a)' at graph 0.785, graph -0.08\n"
                        "set title 'Kernel Splitting for %d-level MSM'\n"
                        "set grid\n"
//                      "set style data lines\n"
                        "set style data linespoints\n"
                        "set yrange [ 0.0 : %7.3f ]\n"
                        "plot ";
    char*   buf2_2 =    "data_file using 1:%d title \"g_%d\" lc rgb \"black\",";
    char*   buf2_3 =    "data_file using 1:%d title \"1/r\" lc rgb \"black\"\n"
                        "pause -1\n"
                        "quit\n";

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
	// Adjust F[0][0] and F[nlev][0] to be large number (Not necessary on OSX)
//	F[0][0] = 1000.0;
//	F[nlev][0] = 1000.0;
	for (i = 0; i <= samples; i++)
	{
// FIXME - There is probably a better way to format this file... binary data would be most accurate, right?
		// Add X to buffer
		buflen = sprintf(buf, "%.32f", X[i]);
		assert(bufmax > buflen);

		for (l = 0; l < nlev; l++)
		{
			// Add F[l][i] to buffer
			buflen += sprintf(&buf[buflen], " %.32f", F[l][i]);
			assert(bufmax > buflen);
		}
		// Add F[nlev][i] to buffer
		buflen += sprintf(&buf[buflen], " %.32f\n", F[nlev][i]);
		assert(bufmax > buflen);

		// Write buffer to file
		bytes = fwrite(buf, sizeof(char), buflen, data);
		if (bytes < buflen)
		{
			printf("<%lu> bytes written to temporary file <%s> instead of <%lu>\n", bytes, data_file, buflen);
			break;
		}

		// Clear out buffer for next usage
		memset(buf, 0, bufmax+1);
	}

	// Close DATA file (ensure buffer is flushed to disk)
	if (fclose(data))
	{
		printf("Error closing DATA file.\n");
	}

	// Create COMMAND file to plot all of the columns above
	cmd = fopen(cmd_file, "w");
	assert(cmd != NULL);

// Write CMD file
	buf2max = strlen(buf2_1) + nlev*strlen(buf2_2) + strlen(buf2_3) + 64;
	buf2 = (char*) dynvec(buf2max+1,sizeof(char));

	// Create COMMAND file buffer
	buflen = sprintf(buf2, buf2_1, GP_TERM, nlev, 1.33*F[1][0]);
	assert(buf2max > buflen);

	// Plot a line for each level (columns 1 -> nlev+1)
	for (i = 0; i < nlev; i++)
	{
		buflen += sprintf(&buf2[buflen], buf2_2, i+2, i);
		assert(buf2max > buflen);
	}

	// Plot line for 1/r (column nlev+2) and finish file
	buflen += sprintf(&buf2[buflen], buf2_3, nlev+2);
	assert(buf2max > buflen);

	// Write buffer to file
	bytes = fwrite(buf2, sizeof(char), buflen, cmd);
	if (bytes < buflen)
	{
		printf("<%lu> bytes written to temporary file <%s> instead of <%lu>\n", bytes, cmd_file, buflen);
	}

	// Close CMD file (ensure buffer is flushed to disk)
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

	// Delete command data file
	if (remove(cmd_file))
	{
		printf("Error removing COMMAND file <%s>.\n", cmd_file);
	}

	// Free dynamically allocated memory
	dynfree(data_file);
	dynfree(cmd_file);
	dynfree(buf);
	dynfree(buf2);
}

/*
Tests gamma, theta, and theta_star without regard to method or scaling
*/
void gamma_test_all(void)
{
	int			i = 0;
	int			samples = 0;
	int			k = 0;
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	char*		data_file = NULL;
	double*		c = NULL;

	/**************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%d", &k);
	assert(k > 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));
	c = (double*) dynvec(k+1,sizeof(double));
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + GAMMA_DATA_LEN + 1,
							   sizeof(char));

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
	// Test Theta* - short range splitting function
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = theta_star(c, k, X[i], &DF[i]);
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
		F[i] = theta(c, k, X[i], &DF[i]);
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

	// Test Gamma - long range splitting function (infinite)
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
}

// End of file