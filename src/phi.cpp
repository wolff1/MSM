//-------|---------|---------|---------|---------|---------|---------|---------|
/*
centered b-spline phi
*/

#include "stdafx.h"

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(int p, double x, double* dphi)
{
	double f = 0.0;
	double df = 0.0;
	double** Q = NULL;
	double** Qp = NULL;
	double u = 0.0;
	short k = 0;
	short v = 0;

	if (fabs(x) < (double)(p/2.0))
	{
		// Dynamically create memory for Q and Qp 2D arrays
		Q = dynarr_d(p,p);		// Q(u)
		Qp = dynarr_d(p,p);		// Q'(u)
		u = x + (double)p/2.0;	// b/c centered B-spline

		// Start with k = 1 (Q_1 is the indicator function on [0,1[)
		k = 1;
		for (v = 0; v <= p-k; v++)
		{
			Q[0][v] = ((u-(double)v) < 1.0 && (u-(double)v) >= 0.0);
		}

		// use recurrance to build up to k=p-1
		for (k = 2; k <= p; k++)
		{
			for (v = 0; v <= p-k; v++)
			{
				Qp[k-1][v] = Q[k-2][v] - Q[k-2][v+1];
				Q [k-1][v] = (k*Q[k-2][v+1] + (u-(double)v)*Qp[k-1][v])/(k-1);
			}
		}

		// Set return variables
		f = Q[p-1][0];
		df = Qp[p-1][0];

		// Free dynamically allocated memory for Q and Qp 2D arrays
		dynfree(Q[0]);		// Free the data buffer memory
		dynfree(Qp[0]);
		dynfree(Q);			// Free the pointer-pointer memory
		dynfree(Qp);
	}

	if (dphi != NULL)
		*dphi = df;
	return f;
}

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of B-spline
g2fg is (p/2 + 1)-vector
*/
void compute_g2fg(short p, double* g2fg)
{
	short		n = 0;
	short		p_2 = (short)p/2;
	long		num = p;			// numerator
	long		den = p_2;			// denominator

	assert(g2fg != NULL);

	// Account for 2^(1-p) in denominator once
	for (n = 1; n < p; n++)
	{
		den *= 2;
	}

	// Build initial state of "p choose p/2"
	for (n = 1; n < p_2; n++) // exclude p_2 b/c it is in num and den
	{
		num *= (p-n);	// [(p)(p-1)(p-2)...(p/2+1)]
		den *= n;		// [(p/2)(p/2-1)...(2)(1)]*[2^(p-1)]
	}
	g2fg[0] = (double) num / (double) den;

	// Compute [p choose p/2 + n] with 1 <= n <= p/2
	for (n = 1; n <= p_2; n++)
	{
		g2fg[n] = (g2fg[n-1]*(double)(p_2-n+1)) / (double)(p_2+n);
	}
}

/*** DRIVER FUNCTIONS BELOW ***/

/*
Tests phi and phi' without regard to method or scaling.
*/
void phi_test_all(void)
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
void print_nesting_coefficients(void)
{
	short		p = 0;
	double*		J = NULL;

	// Get p where p-1 is degree of interpolant
	printf("Please enter accuracy parameter p = ");
	scanf("%hd", &p);
	assert(p % 2 == 0);

	J = (double*) dynvec(p/2+1,sizeof(double));
	assert(J != NULL);

	// Calculate the coefficients to display
	compute_g2fg(p, J);

	for (short i = 0; i <= p/2; i++)
	{
		printf("J[%hd] = %lf\n", i, J[i]);
	}

	// Free dynamically allocated memory
	dynfree(J);
}

/*
*/
void phi_nesting_test(void)
{
	// Static memory variables
	int			i = 0;
	int			j = 0;
	int			samples = 0;
	int			p = 0;
	char*		cmd_file_name = "nesting_";
	char*		cmd_file_extension = ".gp";
	size_t		cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN + strlen(cmd_file_extension);
	size_t		bufmax = 0;
	size_t		buflen = 0;
	size_t		bytes = 0;
	double		x = 0.0;
	double		ff = 0.0;	// Function value on fine grid
	double		fc = 0.0;	// Function value on course grid
	double		f = 0.0;	// Sum of fine grid contributions
	double		maxf = 0.0;
	double		tol = pow(0.1,15);
    char*		buf2_1 =    "set term %s\n"
							"set termoption dash\n"
		                    "set xlabel 'u'\n"
//			                "set ylabel 'g_l(r)'\n" // below two lines rotates by -90 degrees
							"set lmargin 10\n"
							"set label 1 'f(u)' at graph -0.1, graph 0.5\n"
							"set title 'Nesting of Spline Function Spaces'\n"
							"set grid\n"
							"set style line 1 lc rgb \"cyan\"\n"
							"set style line 2 lc rgb \"blue\"\n"
							"set style line 3 lc rgb \"red\"\n"
//							"set style data lines\n"
//							"set style data linespoints\n"
							"set yrange [ 0.0 : %7.3f ]\n"
							"plot ";
    //char*		buf2_2 =    "data_file using 1:%d title \"Phi(2u%+d)\" lc rgb \"black\",";
    //char*		buf2_3 =    "data_file using 1:%d title \"Sum of fine grids\" lc rgb \"black\","
				//			"data_file using 1:%d title \"Phi(u)\" lc rgb \"black\"\n"
				//			"pause -1\n"
				//			"quit\n";
    char*		buf2_2 =    "data_file using 1:%d with lines title \"Phi(2u%+d)\" ls 1,";
    char*		buf2_3 =    "data_file using 1:%d with lines title \"Sum of fine grids\" ls 2,"
							"data_file using 1:%d with lines title \"Phi(u)\" ls 3\n"
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
	g2fg = (double*) dynvec(p/2+1,sizeof(double));
	assert(data_file != NULL);
	assert(cmd_file != NULL);
	assert(g2fg != NULL);

	// Calculate the coefficients to display
	compute_g2fg(p, g2fg);

	// Build file name(s)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(cmd_file, "%s%s%s%s", GP_CMD_DIR, cmd_file_name, GP_TERM, cmd_file_extension);

	//	Create DATA file
	fp = fopen(data_file, "w");
	assert(fp != NULL);

	bufmax = 64*(p+4);	// p+4 columns, each a max of 64 chars wide 
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
		for (j = 1; j <= p/2; j++)
		{
			ff = g2fg[j]*phi(p, 2.0*x-j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

		ff = g2fg[0]*phi(p, 2.0*x, NULL);
f += ff;
		buflen += sprintf(&buf[buflen], " %.32f", ff);
		assert(bufmax > buflen);

		for (j = 1; j <= p/2; j++)
		{
			ff = g2fg[j]*phi(p, 2.0*x+j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

		// Column p+3 (sum of fine grid interpolant values)
		buflen += sprintf(&buf[buflen], " %.32f", f);
		assert(bufmax > buflen);

		// Column p+4 (course grid interpolant values)
		fc = phi(p, x, NULL);
		buflen += sprintf(&buf[buflen], " %.32f\n", fc);
		assert(bufmax > buflen);
		if (fc > maxf)
		{
			maxf = fc;
		}
if (fabs(f - fc) >= tol)
{
	printf("|f-fc| = %e\n", fabs(f-fc));
}
assert(fabs(f - fc) < tol);

		// Write buffer to file
		bytes = fwrite(buf, sizeof(char), buflen, fp);
		if (bytes < buflen)
		{
			printf("<%lu> bytes written to temporary file <%s> instead of <%lu>\n", bytes, data_file, buflen);
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
	bufmax = strlen(buf2_1) + (p+2)*(strlen(buf2_2)+10) + strlen(buf2_3) + 64;
	buf = (char*) dynvec(bufmax+1,sizeof(char));

	/*****WRITE THE COMMAND FILE TO DISPLAY THE RESULTS OF ABOVE *************/
	fp = fopen(cmd_file, "w");
	assert(fp != NULL);

	// Create COMMAND file buffer
	buflen = sprintf(buf, buf2_1, GP_TERM, 1.25*maxf);
	assert(bufmax > buflen);

	// Plot a line for each fine grid b-spline (p+1 fine grids, 1 sum of fine grids)
	for (i = 0; i < p+1; i++)
	{
		buflen += sprintf(&buf[buflen], buf2_2, i+2, i-p/2);
		assert(bufmax > buflen);
	}

	// Plot line for 1/r (column nlev+2) and finish file
	buflen += sprintf(&buf[buflen], buf2_3, p+3, p+4);
	assert(bufmax > buflen);

	// Write buffer to file
	bytes = fwrite(buf, sizeof(char), buflen, fp);
	if (bytes < buflen)
	{
		printf("<%lu> bytes written to temporary file <%s> instead of <%lu>\n", bytes, cmd_file, buflen);
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

//	End of file