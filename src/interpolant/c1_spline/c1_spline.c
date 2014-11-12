//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.c - C1 phi (CUBIC ONLY)
*/

#include "all.h"

// These are for the CUBIC "shifted-powers" version
double c[2][4] = {{0.0, -0.5, 2.0, 1.5},
				  {0.0, 0.0, -0.5,-0.5}};

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phiC1(int p, double x, double* dphi)
{
	double	f = 0.0;
	double	df = 0.0;
	double	signx = 0.0;
	short	piece = 0;
	short	i = 0;

	if (x > 0.0)
		signx = 1.0;
	else if (x < 0.0)
		signx = -1.0;

	x = fabs(x);
	if (x < p/2)
	{
		piece = (short)floor(x);
		x = x - (piece + 1.0);
		f = c[piece][p-1];
		df = c[piece][p-1]*(p-1);
		for (i = p-2; i > 0; i--)
		{
			f = f*x + c[piece][i];
			df = df*x + c[piece][i]*(i);
		}
		f = f*x + c[piece][0];
		df = df*signx;
	}

	if (dphi != NULL)
		*dphi = df;
	return f;
}

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of interpolant
g2fg is (p + 1)-vector
*/
void compute_g2fgC1(short p, double* g2fg)
{
	short		n = 0;

	assert(g2fg != NULL);

	for (n = 0; n <= p; n++)
	{
		g2fg[n] = phiC1(p,n/2.0,NULL);
	}
}

/*** DRIVER FUNCTIONS BELOW ***/

/*
Tests phi and phi' without regard to method or scaling.
*/
void phi_test_allC1(void)
{
	int			i = 0;
	int			samples = 0;
	int			p = 0;
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;

	/**************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%d", &p);
	assert(p % 2 == 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));
	assert(X != NULL);
	assert(F != NULL);
	assert(DF != NULL);

	/**************************************************************************/
	// Test interpolant where p must be even
	for (i = 0; i <= samples; i++)
	{
		X[i] = (double)(-p/2.0) + ((double)i/(double)samples)*((double)p);
		F[i] = phiC1(p, X[i], &DF[i]);
/*
		printf("%02d:\tphi(%f) = %f\tphi'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Phi(x)
	printf("Plotting Phi(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, F);
	//	Show Phi'(x)
	printf("Plotting Phi'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	/**************************************************************************/
	// Free allocated memory
	dynfree(X);
	dynfree(F);
	dynfree(DF);
}

/*
Simple interface to output nesting coefficients to user.
*/
void print_nesting_coefficientsC1(void)
{
	short		i = 0;
	short		p = 0;
	double*		J = NULL;

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%hd", &p);
	assert(p % 2 == 0);

	J = (double*) dynvec(p+1,sizeof(double));
	assert(J != NULL);

	// Calculate the coefficients to display
	compute_g2fgC1(p, J);

	for (i = 0; i <= p; i++)
	{
		printf("J[%hd] = %lf\n", i, J[i]);
	}

	// Free dynamically allocated memory
	dynfree(J);
}

/*
*/
void phi_nesting_testC1(void)
{
	// Static memory variables
	int			i = 0;
	int			j = 0;
	int			samples = 0;
	int			p = 0;
	char*		cmd_file_name = "nesting_";
	char*		cmd_file_extension = ".gp";
	size_t		cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN +
									strlen(cmd_file_extension);
	size_t		bufmax = 0;
	size_t		buflen = 0;
	size_t		bytes = 0;
	double		x = 0.0;
	double		ff = 0.0;	// Function value on fine grid
	double		fc = 0.0;	// Function value on course grid
	double		f = 0.0;	// Sum of fine grid contributions
	double		maxf = 0.0;
	double		minf = 0.0;
	double		tol = pow(0.1,15);
    char*		buf2_1 =    "set term %s\n"
							"set termoption dash\n"
		                    "set xlabel 'u'\n"
			                "set ylabel 'f(u)'\n"
							// below two lines rotates by -90 degrees
							//"set lmargin 10\n"
							//"set label 1 'f(u)' at graph -0.1, graph 0.5\n"
							"set title 'Non-nesting of C1 Function Spaces'\n"
							"set grid\n"
							"set style data lines\n"
							"set yrange [ %7.3f : %7.3f ]\n"
							"set key box\n"
							"plot ";
    char*		buf2_2 =    "data_file using 1:%d with lines title "
							"\"Phi(2u%+d)\",";
    char*		buf2_3 =    "data_file using 1:%d with lines title "
							"\"Sum of fine grids\","
							"data_file using 1:%d with lines title "
							"\"Phi(u)\"\n"
							"pause -1\n"
							"quit\n";
	// Dynamically allocated memory variables
	FILE*		fp = NULL;
	char*		data_file = NULL;
	char*		cmd_file = NULL;
	char*		buf = NULL;
	double*		g2fg = NULL;

	/**************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%d", &p);
	assert(p % 2 == 0);

	// Create arrays for dependent and independent variables
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + GP_DATA_TMP_LEN + 1,sizeof(char));
	cmd_file = (char*) dynvec(GP_CMD_DIR_LEN + cmd_file_name_len + 1,sizeof(char));
	g2fg = (double*) dynvec(p+1,sizeof(double));
	assert(data_file != NULL);
	assert(cmd_file != NULL);
	assert(g2fg != NULL);

	// Calculate the coefficients to display
	compute_g2fgC1(p, g2fg);

	// Build file name(s)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(cmd_file, "%s%s%s%s",
					GP_CMD_DIR, cmd_file_name, GP_TERM, cmd_file_extension);

	//	Create DATA file
	fp = fopen(data_file, "w");
	assert(fp != NULL);

	bufmax = 64*(p+4);	// 2p+4 columns, each a max of 64 chars wide 
	buf = (char*) dynvec(bufmax+1,sizeof(char));	// + 1 for NULL

	/**** WRITE THE DATA FILE WHILE CALCULATING ALL THE REQUIRED VALUES*******/
	for (i = 0; i <= samples; i++)
	{
		//	Column 1 (Independent variable)
		x = (double)(-p/2.0) + ((double)i/(double)samples)*((double)p);

		buflen = sprintf(buf, "%.32f", x);
		assert(bufmax > buflen);
f = 0.0;

		//	Columns 2 to p+2 (fine grid interpolant values)
		for (j = p-1; j > 0; j-=2)
		{
			ff = g2fg[j]*phiC1(p, 2.0*x-j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

		ff = g2fg[0]*phiC1(p, 2.0*x, NULL);
f += ff;
		buflen += sprintf(&buf[buflen], " %.32f", ff);
		assert(bufmax > buflen);

		for (j = 1; j < p; j+=2)
		{
			ff = g2fg[j]*phiC1(p, 2.0*x+j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

		// Column p+3 (sum of fine grid interpolant values)
		buflen += sprintf(&buf[buflen], " %.32f", f);
		assert(bufmax > buflen);

		// Column p+4 (course grid interpolant values)
		fc = phiC1(p, x, NULL);
		buflen += sprintf(&buf[buflen], " %.32f\n", fc);
		assert(bufmax > buflen);
		if (fc > maxf)
		{
			maxf = fc;
		}
		if (fc < minf)
		{
			minf = fc;
		}
/*
if (fabs(f - fc) >= tol)
{
	printf("|f-fc| = %e\n", fabs(f-fc));
}
assert(fabs(f - fc) < tol);
*/
		// Write buffer to file
		bytes = fwrite(buf, sizeof(char), buflen, fp);
		if (bytes < buflen)
		{
			printf("<%lu> bytes written to temp file <%s> instead of <%lu>\n",
					bytes, data_file, buflen);
			break;
		}

		// Clear out buffer for next usage
		memset(buf, 0, bufmax+1);
	}

	// Close DATA file (ensure buffer is flushed to disk)
	if (fclose(fp))
	{
		printf("Error closing DATA file.\n");
	}

	// Free buffer and re-allocate it for next file
	dynfree(buf);
	bufmax = strlen(buf2_1) + (2*p+2)*(strlen(buf2_2)+10) + strlen(buf2_3) + 64;
	buf = (char*) dynvec(bufmax+1,sizeof(char));

	/*****WRITE THE COMMAND FILE TO DISPLAY THE RESULTS OF ABOVE *************/
	fp = fopen(cmd_file, "w");
	assert(fp != NULL);

	// Create COMMAND file buffer
	buflen = sprintf(buf, buf2_1, GP_TERM, -0.2, 1.2);
	assert(bufmax > buflen);

	// Plot a line for each fine grid (p+1 fine grids, 1 sum of fine grids)
	for (i = 0; i < p/2; i++)
	{
		buflen += sprintf(&buf[buflen], buf2_2, i+2, 2*i-3);
		assert(bufmax > buflen);
	}
	buflen += sprintf(&buf[buflen], buf2_2, p/2+2, 0);
	assert(bufmax > buflen);
	for (i = 0; i < p/2; i++)
	{
		buflen += sprintf(&buf[buflen], buf2_2, p/2+3+i, 2*i+1);
		assert(bufmax > buflen);
	}

	// Plot line for 1/r (column nlev+2) and finish file
	buflen += sprintf(&buf[buflen], buf2_3, p+3, p+4);
	assert(bufmax > buflen);

	// Write buffer to file
	bytes = fwrite(buf, sizeof(char), buflen, fp);
	if (bytes < buflen)
	{
		printf("<%lu> bytes written to temp file <%s> instead of <%lu>\n",
				bytes, cmd_file, buflen);
	}

	// Close COMMAND file (ensure buffer is flushed to disk)
	if (fclose(fp))
	{
		printf("Error closing COMMAND file.\n");
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

	/**************************************************************************/
	// Free allocated memory
	dynfree(data_file);
	dynfree(cmd_file);
	dynfree(buf);
	dynfree(g2fg);
}

/*
*/
void driverC1(void)
{
//	phi_test_allC1();
//	print_nesting_coefficientsC1();
	phi_nesting_testC1();
}

//	End of file