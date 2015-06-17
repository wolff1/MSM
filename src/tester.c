//-------|---------|---------|---------|---------|---------|---------|---------|
/*
tester.c - 
*/

//	NOTE:		CURRENTLY UNUSED!

#include "all.h"
#include "tester.h"

#if 0
//	**** POLYNOMIAL ****
void test_mpoly(void)
{
	double*		a = NULL;
	double*		b = NULL;
	double*		c = NULL;
	short		i = 0;
	short		d1 = 0;
	short		d2 = 0;
	short		d3 = 0;

	//	Get degrees from user
	printf("Please enter degree 1: ");
	scanf("%hd", &d1);
	printf("Please enter degree 2: ");
	scanf("%hd", &d2);
	d3 = d1+d2;

	// Dynamically allocate memory
	a = (double*) dynvec(d1+1, sizeof(double));
	b = (double*) dynvec(d2+1, sizeof(double));
	c = (double*) dynvec(d3+1, sizeof(double));

	// Initialize a
	for (i = 0; i <= d1; i++)
	{
		a[i] = (double)i + 1.0;
	}

	// Initialize b
	for (i = 0; i <= d2; i++)
	{
		b[i] = (double)-i - 1.0;
	}

	mpoly(d1, a, d2, b, c);

	printf("\n");
	for (i = 0; i <= d3; i++)
	{
		printf("%f ", c[i]);
	}
	printf("\n");

	// Free dynamically allocated memory
	dynfree(a);
	dynfree(b);
	dynfree(c);
}

//	**** INTERPOLANT ****
void phi_test_all(void)
{
#if 0
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
#endif
    //  Compute fixed Phi coefficients for polynomial evaluation
    short       p = 0;
    short       p_2 = 0;
    double**    g2p = NULL;
    double**    A = NULL;
    double*     b = NULL;
    short       i = 0;
    short       j = 0;
    short       k = 0;
    double      x = 0.0;
    double      xmk = 0.0;
    double      dx = 0;
	lapack_int  rc = 0;
	lapack_int* piv = NULL;

    short       samples = 0;
    double*     y1 = NULL;
    double*     dy1 = NULL;
    double*     y2 = NULL;
    double*     dy2 = NULL;
    double*     xs = NULL;

	double		max = 0.0;
	double		min = 1.0;
	double		maxd = 0.0;
	double		mind = 1.0;

    printf("Enter p: ");
    scanf("%hd", &p);

    p_2 = p/2;
    dx = 1.0/(p_2-2+1);
    g2p = (double**) dynarr_d(p_2,p);
    A = (double**) dynarr_d(p,p);
    b = (double*) dynvec(p, sizeof(double));

    for (k = 0; k < p_2; k++)
    {
        x = k;
        for (j = 0; j < p_2; j++)
        {
            //  down column -> different values of x
//            printf("k=%hd, j=%hd, x=%f\n", k, j, x);
            xmk = x - (k+1);

            //  across powers in row
            A[j][0] = 1.0;
            A[j+p_2][0] = 0.0;
            A[j+p_2][1] = 1.0;
            for (i = 1; i < p-1; i++)
            {
                A[j][i] = A[j][i-1]*xmk;
                A[j+p_2][i+1] = (i+1)*A[j][i];
            }
            A[j][p-1] = A[j][p-2]*xmk;

            //  right hand side of linear system
            b[j] = phi(p, x, &b[j+p_2]);
            x += dx;
        }

//        display_dynarr_d(A, p, p);
//        display_vector_d(b, p);

        //  Solve Ax = b with x being g2p[k]
        // NOTE: b is overwritten with x, A is overwritten with LU
        piv = (lapack_int*) dynvec(p,sizeof(lapack_int));
        rc = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) p, (lapack_int) 1,
                            A[0], (lapack_int) p, piv, b, (lapack_int) 1);
        dynfree(piv);
        assert(rc == 0); // zero is SUCCESS

        for (i = 0; i < p; i++)
        {
            g2p[k][i] = b[i];
        }

//        display_vector_d(g2p[k], p);
    }

    //  Test new phi
    printf("Enter samples: ");
    scanf("%hd", &samples);

    y1 = (double*) dynvec(samples+1, sizeof(double));
    dy1 = (double*) dynvec(samples+1, sizeof(double));
    y2 = (double*) dynvec(samples+1, sizeof(double));
    dy2 = (double*) dynvec(samples+1, sizeof(double));
    xs = (double*) dynvec(samples+1, sizeof(double));

	for (i = 0; i <= samples; i++)
    {
        xs[i] = (double)i*p/samples - (double)p_2;
        y1[i] = phi(p,xs[i],&dy1[i]);
    }

    new_phi(p, g2p, samples+1, xs, y2, dy2);

	//	Determine min/max (absolute) errors
    for (i = 0; i <= samples; i++)
    {
		if (fabs(y2[i]-y1[i]) > max)
			max = fabs(y2[i]-y1[i]);
		if (fabs(y2[i]-y1[i]) < min)
			min = fabs(y2[i]-y1[i]);

		if (fabs(dy2[i]-dy1[i]) > maxd)
			maxd = fabs(dy2[i]-dy1[i]);
		if (fabs(dy2[i]-dy1[i]) < mind)
			mind = fabs(dy2[i]-dy1[i]);
	}
	printf("min/max y = (%e,%e), dy= (%e,%e)\n", min,max, mind,maxd);
    
    dynfree(g2p[0]);
    dynfree(g2p);
    dynfree(A[0]);
    dynfree(A);
    dynfree(b);
	dynfree(y1);
	dynfree(dy1);
	dynfree(y2);
	dynfree(dy2);
	dynfree(xs);
}

void print_nesting_coefficients(void)
{
	short		i = 0;
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

	for (i = 0; i <= p/2; i++)
	{
		printf("J[%hd] = %lf\n", i, J[i]);
	}

	// Free dynamically allocated memory
	dynfree(J);
}
#endif

void phi_nesting_test(void)
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
							"set title 'Nesting of Spline Function Spaces'\n"
							"set grid\n"
							//"set style line 1 lc rgb \"cyan\"\n"
							//"set style line 2 lc rgb \"blue\"\n"
							//"set style line 3 lc rgb \"red\"\n"
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
	g2fg = (double*) dynvec(p/2+1,sizeof(double));
	assert(data_file != NULL);
	assert(cmd_file != NULL);
	assert(g2fg != NULL);

	// Calculate the coefficients to display
//	compute_g2fg(p, g2fg);

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
		for (j = p/2; j >= 1; j--)
		{
//			ff = g2fg[j]*phi(p, 2.0*x-j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

//		ff = g2fg[0]*phi(p, 2.0*x, NULL);
f += ff;
		buflen += sprintf(&buf[buflen], " %.32f", ff);
		assert(bufmax > buflen);

		for (j = 1; j <= p/2; j++)
		{
//			ff = g2fg[j]*phi(p, 2.0*x+j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

		// Column p+3 (sum of fine grid interpolant values)
		buflen += sprintf(&buf[buflen], " %.32f", f);
		assert(bufmax > buflen);

		// Column p+4 (course grid interpolant values)
//		fc = phi(p, x, NULL);
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
	buflen = sprintf(buf, buf2_1, GP_TERM, 1.25*minf, 1.25*maxf);
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
#if 0
//	**** B_SPLINE ****
void test_blurring_operator(void)
{
	short	p = 0;
	short	p_2 = 0;
	short	k = 0;
	double*	B = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);
	p_2 = p/2;

	// Dynamically allocate memory for blurring operator
	B = (double*) dynvec(p_2,sizeof(double));

	b_spline_compute_blurring_operator(p_2-1, B);

	// Print c_n
	printf("B = (%f", B[0]);
	for (k = 1; k < p_2; k++)
	{
		printf(", %f", B[k]);
	}
	printf(")\n");

	// Free dynamically allocated memory
	dynfree(B);
}

void test_omega(void)
{
	short		i = 0;
	short		p = 0;
	short		p_2 = 0;
	short		n = 0;
	double*		omega = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);
	p_2 = p/2;
	printf("Please enter n: ");
	scanf("%hd", &n);

	// Dynamically allocate memory
	omega = (double*) dynvec(n, sizeof(double));

	compute_omega(p, n, omega);

	printf("\n");
	for (i = 0; i < n; i++)
	{
		printf("%02hd: %f\n", i, omega[i]);
	}
	printf("\n");

	// Free dynamically allocated memory
	dynfree(omega);
}

void test_omegap(void)
{
	short		i = 0;
	short		p = 0;
	short		p_2 = 0;
	short		mu = 0;
	double*		omegap = NULL;

	printf("Please enter p: ");
	scanf("%hd", &p);
	p_2 = p/2;
	printf("Please enter mu: ");
	scanf("%hd", &mu);

	// Dynamically allocate memory
	omegap = (double*) dynvec(p_2+mu+1, sizeof(double));

	compute_omega_prime(p, mu, omegap);

	printf("\n");
	for (i = 0; i <= p_2+mu; i++)
	{
		printf("%02hd: %f\n", i, omegap[i]);
	}
	printf("\n");

	// Free dynamically allocated memory
	dynfree(omegap);
}

void test_convert_to_shifts(void)
{
	double*		cd = NULL;
	double*		s = NULL;
	short		degree = 2;
	short		i = 0;
	short		j = 0;

	for (i = 0; i <= 4; i++)
	{
		// Allocate
		cd = (double*) dynvec(i+1, sizeof(double));
		s  = (double*) dynvec(i+1, sizeof(double));

		// set one term and display
		cd[i] = 1.0;
		for (j = 0; j <= i; j++)
		{
			printf("cd[%02hd] = %f\n", j, cd[j]);
		}
		printf("\n");

		// Convert cd to s
		b_spline_convert_delta2_to_shifts(i, cd, s);

		// display result
		for (j = 0; j <= i; j++)
		{
			printf("s[%02hd] = %f\n", j, s[j]);
		}
		printf("\n----------------------------------\n");

		// Free
		dynfree(cd);
		dynfree(s);
	}
}

//	**** C1 SPLINE ****
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
#endif
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
//	compute_g2fgC1(p, g2fg);

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
//			ff = g2fg[j]*phiC1(p, 2.0*x-j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

//		ff = g2fg[0]*phiC1(p, 2.0*x, NULL);
f += ff;
		buflen += sprintf(&buf[buflen], " %.32f", ff);
		assert(bufmax > buflen);

		for (j = 1; j < p; j+=2)
		{
//			ff = g2fg[j]*phiC1(p, 2.0*x+j, NULL);
f += ff;
			buflen += sprintf(&buf[buflen], " %.32f", ff);
			assert(bufmax > buflen);
		}

		// Column p+3 (sum of fine grid interpolant values)
		buflen += sprintf(&buf[buflen], " %.32f", f);
		assert(bufmax > buflen);

		// Column p+4 (course grid interpolant values)
//		fc = phiC1(p, x, NULL);
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

void driverC1(void)
{
//	phi_test_allC1();
//	print_nesting_coefficientsC1();
	phi_nesting_testC1();
}

//	**** EVEN_POWERS ****
void gamma_test_all(void)
{
	int			i = 0;
	int			samples = 0;
	short		k = 0;
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	double*		c = NULL;

	size_t		Size = 0;
	void*		Init = NULL;
	void*		Ptr = NULL;

	/**************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%hd", &k);
	assert(k > 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));
	c = (double*) dynvec(k+1,sizeof(double));
	assert(X != NULL);
	assert(F != NULL);
	assert(DF != NULL);
	assert(c != NULL);

	//	Initialize SOFTENER
	if (1)
	{	//	EVEN_POWERS softening (aka "Taylor")
		Size = sizeof(EVEN_POWERS);
		Init = &even_powers_initialize;
	}
	Ptr = (SOFTENER*) dynmem(Size);
	softener_initialize(Ptr, Init, k);

	/**************************************************************************/
	// Test softener
//	gamma_init(k, c);
	for (i = 0; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
//		F[i] = _gamma(c, k, X[i], &DF[i]);
/*
		printf("%02d:\tgamma(%f) = %f\tgamma'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}

	((SOFTENER*)Ptr)->soften(Ptr, samples+1, X, F, DF);

	//	Show Gamma(x)
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, F);
	//	Show Gamma'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);
#if 0
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
	printf("Plotting theta*(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &F[1]);
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
	printf("Plotting theta(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &F[1]);
	//	Show Theta'(x)
	printf("Plotting theta'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	// Test Gamma - long range splitting function (infinite)
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = _gamma(c, k, X[i], &DF[i]);
/*
		printf("%02d:\tgamma(%f) = %f\tgamma'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show gamma(x)
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &F[1]);
	//	Show gamma'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);
#endif
	/**************************************************************************/
	// Free allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F);
	dynfree(DF);
}
#if 0
void plot_splittings(int samples, int nlev, double a, double d,
						double* X, double** F)
{
	// Static memory variables
	char*	cmd_file_name = "splitting_";
	char*	cmd_file_extension = ".gp";
	size_t	cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN +
								strlen(cmd_file_extension);
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
						"set ylabel 'g_l(r)'\n"
						// below two lines rotates by -90 degrees
                        //"set lmargin 10\n"
                        //"set label 1 'g_l(r)' at graph -0.1, graph 0.5\n"
                        "set label 2 '(a)' at graph 0.39, graph -0.08\n"
                        "set label 3 '(2a)' at graph 0.785, graph -0.08\n"
                        "set title 'Kernel Splitting for %d-level MSM'\n"
                        "set grid\n"
//                      "set style data lines\n"
                        "set style data linespoints\n"
                        "set yrange [ 0.0 : %7.3f ]\n"
						"set key box\n"
                        "plot ";
    char*   buf2_2 =    "data_file using 1:%d title \"g_%d\" lc rgb \"black\",";
    char*   buf2_3 =    "data_file using 1:%d title \"1/r\" lc rgb \"black\"\n"
                        "pause -1\n"
                        "quit\n";

	assert(X != NULL);
	assert(F != NULL);

	// Dynamically allocate memory
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + GP_DATA_TMP_LEN + 1,
								sizeof(char));
	cmd_file = (char*) dynvec(GP_CMD_DIR_LEN + cmd_file_name_len + 1,
								sizeof(char));

	// Build file name(s)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(cmd_file, "%s%s%s%s",
					GP_CMD_DIR, cmd_file_name, GP_TERM, cmd_file_extension);

	// Create DATA file with X, F[0], F[1], ..., F[nlev] columns and samples rows
	data = fopen(data_file, "w");
	assert(data != NULL);

	// Write DATA file
	bufmax = 64*(nlev+2);	// nlev + 2 columns, each a max of 64 chars wide 
	buf = (char*) dynvec(bufmax+1,sizeof(char));	// + 1 for NULL
	// Adjust F[0][0] and F[nlev][0] to be large number (Not necessary on OSX)
	F[0][0] = 1000.0;
	F[nlev][0] = 1000.0;
	for (i = 0; i <= samples; i++)
	{
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
			printf("<%lu> bytes written to temp file <%s> instead of <%lu>\n",
				bytes, data_file, buflen);
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
		printf("<%lu> bytes written to temporary file <%s> instead of <%lu>\n",
				bytes, cmd_file, buflen);
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
#endif
void splitting_test(void)
{
	// Static memory variables
	int			i = 0;
	int			samples = 0;
	short		k = 0;
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

	// Get nlev where nlev is the number of levels
	printf("Please enter levels parameter nlev = ");
	scanf("%d", &nlev);
	assert(nlev > 1);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%hd", &k);
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
	data_file = (char*) dynvec(GP_DATA_DIR_LEN +
						MAX(PHI_DATA_LEN,GAMMA_DATA_LEN) + 1, sizeof(char));
	assert(X != NULL);
	assert(F != NULL);
	assert(DF != NULL);
	assert(c != NULL);
	assert(data_file != NULL);

	// Initialize gamma coefficients
//	gamma_init(k, c);

	// Evaluate different splitting functions over domain
	for (i = 0; i <= samples; i++)
	{
f = 0.0;	// sanity check
df = 0.0;	// sanity check
		X[i] = (d*(double)i/(double)samples);

		// Theta* - Short range part of splitting (finite)
		al = one_over_a;
//		F[0][i] = al*theta_star(c, k, al*X[i], &DF[0][i]);
		DF[0][i] = al*al*DF[0][i];
f += F[0][i];
df += DF[0][i];

		for (l = 1; l < nlev - 1; l++)
		{
			// Theta - Intermediate long range part(s) of splitting (finite)
//			F[l][i] = al*theta(c, k, al*X[i], &DF[l][i]);
			DF[l][i] = al*al*DF[l][i];
f += F[l][i];
df += DF[l][i];

			al = 0.5*al;
		}
		// Gamma - Top level long range part of splitting (infinit)
//		F[l][i] = _gamma(c, k, al*X[i], &DF[l][i]);
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
		printf("i = %d, X = %f, F = %f, f = %f, |.| = %e\n",
				i, X[i], F[nlev][i], f, Frel);
	if (DFrel >= tol)
		printf("i = %d, X = %f, DF = %f, df = %f, |.| = %e\n",
				i, X[i], DF[nlev][i], df, DFrel);
	assert(Frel < tol);
	assert(DFrel < tol);
}
	}

	// Plot splittings on single graph along with f(x) = 1/x
//	plot_splittings(samples, nlev, a, d, X, F);
////	plot_splittings(samples, nlev, k, a, d, X, DF);

	// Free dynamically allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F[0]);
	dynfree(F);
	dynfree(DF[0]);
	dynfree(DF);
	dynfree(data_file);
}
#if 0
void test_thetas(void)
{
	int		i = 0;
	int		samples = 1000;
	double	t1 = 0.0;
	double	dt1 = 0.0;
	double	t2 = 0.0;
	double	dt2 = 0.0;
	double	e = 0.0;
	double	de = 0.0;
	short	k = 0;
	double*	c = NULL;
	double	x = 0.0;

	printf("Enter samples: ");
	scanf("%d", &samples);

	printf("Enter k: ");
	scanf("%hd", &k);

	c = (double*) dynvec(k+1,sizeof(double));

	// Get gamma coefficients
	gamma_init(k, c);

	for (i = 0; i <= samples; i++)
	{
		x = ((double)i/(double)samples)*2.5;
		t1 = theta(c,k,x,&dt1);
		t2 = thetap(c,k,x,&dt2);

		e = fabs(t1 - t2);
		de = fabs(dt1 - dt2);

		printf("%f: Error = %e, dError = %e\n", x, e, de);
	}

	dynfree(c);
}
#endif

//	**** OTHER ****
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
	char*		cmd_file_name = "sinc_";
	char*		cmd_file_extension = ".gp";
	size_t		cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN +
									strlen(cmd_file_extension);
	FILE*		fp = NULL;
	char*		data_file = NULL;
	char*		cmd_file = NULL;
	char*		buf = NULL;
	size_t		bufmax = 512;//6*32;
	size_t		buflen = 0;
	size_t		bytes = 0;

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
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + GP_DATA_TMP_LEN + 1,sizeof(char));
	cmd_file = (char*) dynvec(GP_CMD_DIR_LEN + cmd_file_name_len + 1,sizeof(char));
	buf = (char*) dynvec(bufmax+1, sizeof(char));

	//	Compute omega values
//	compute_omega(4, omega_max, omega3);
//	compute_omega(6, omega_max, omega5);
//	compute_omega(8, omega_max, omega7);

	// Build file name(s)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(cmd_file, "%s%s%s%s", GP_CMD_DIR, cmd_file_name, GP_TERM, cmd_file_extension);

	//	Create DATA file
	fp = fopen(data_file, "w");
	assert(fp != NULL);

	//	Compute values of interpolants and sinc function
	for (i = 0; i <= samples; i++)
	{
		u = umin + (umax-umin)*(double)i /samples;

		//	Cubic
		p = 4;
		for (j = 0; j < p; j++)
		{
			n = (short)floor(u) - p/2 + 1 + j;
			assert(abs(n) < omega_max);
//			results[0][i] += omega3[abs(n)]*phi(p,u-(double)n,NULL);
		}

		//	Quintic
		p = 6;
		for (j = 0; j < p; j++)
		{
			n = (short)floor(u) - p/2 + 1 + j;
			assert(abs(n) < omega_max);
//			results[1][i] += omega5[abs(n)]*phi(p,u-(double)n,NULL);
		}

		//	Septic
		p = 8;
		for (j = 0; j < p; j++)
		{
			n = (short)floor(u) - p/2 + 1 + j;
			assert(abs(n) < omega_max);
//			results[2][i] += omega7[abs(n)]*phi(p,u-(double)n,NULL);
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

		//	Write data to buffer
		buflen = sprintf(buf, "%.32f %.32f %.32f %.32f %.32f\n", u, results[0][i], results[1][i], results[2][i], results[3][i]);
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

	/*****WRITE THE COMMAND FILE TO DISPLAY THE RESULTS OF ABOVE *************/
	fp = fopen(cmd_file, "w");
	assert(fp != NULL);

	// Write buffer to file (buf set to size of 512 bytes)
	buflen = sprintf(buf,	"set term %s\n"
							"set xlabel 'u'\n"
							"set ylabel 'f(u)'\n"
							"set title 'Fundamental splines of varying degree and sinc function'\n"
							"set style data lines\n"
							"set key box\n"
							"plot [][-0.4:1.2] 0.0 title \"\", "
									"data_file using 1:2 title \"Cubic\", "
										 "\"\" using 1:3 title \"Quintic\", "
										 "\"\" using 1:4 title \"Septic\", "
										 "\"\" using 1:5 title \"Sinc\"\n"
							"pause -1\n"
							"quit\n",
							GP_TERM);
	assert(buflen < bufmax);
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

	//	Free dynamically allocated memory
	dynfree(omega3);
	dynfree(omega5);
	dynfree(omega7);
	dynfree(results[0]);
	dynfree(results);
	dynfree(data_file);
	dynfree(cmd_file);
	dynfree(buf);
}
#if 0
void test_preprocessing(void)
{
	short			k = 0;
	double			a = 0;
	double			h = 0.0;
	double			alpha = 0.0;
	short			p = 0;
	short			p_2 = 0;
	short			mu = 0;
	double*			c = NULL;
	double*			omegap = NULL;
	STENCIL			theta;
	STENCIL			gamma;
	STENCIL			g2g;
	STENCIL			tg2g;

	void*			Init = NULL;
	SOFTENER*		s = NULL;
	
	//	Get parameters from user
	printf("Please enter p: ");
	scanf("%hd", &p);

	printf("Please enter mu: ");
	scanf("%hd", &mu);

	printf("Enter k: ");
	scanf("%hd", &k);

	printf("Enter a: "); // alpha = a/h, radius = 2*alpha
	scanf("%lf", &a);

	printf("Enter h: "); // alpha = a/h, radius = 2*alpha
	scanf("%lf", &h);
/*
	Init = &even_powers_initialize;
	softener_initialize((void**)&s, sizeof(EVEN_POWERS), Init, k);

	//	Compute computed values
	p_2 = p/2;
	alpha = a/h;

	// Dynamically allocate memory
	c = (double*) dynvec(k+1,sizeof(double));
	omegap = (double*) dynvec(p_2+mu+1, sizeof(double));

	//	Pre-pre-processing
	gamma_init(k, c);
	compute_omega_prime(p, mu, omegap);

	//	Pre-processing (Intermediate levels)
	stencil_initialize(&theta, (long) ceil(2.0*alpha), STENCIL_SHAPE_SPHERE);
	stencil_populate(&theta, c, k, STENCIL_FUNCTION_TYPE_THETA, h/a);
	//stencil_display(&theta, h/a);

	stencil_initialize(&g2g, theta.size, theta.shape);
	stencil_shift(&theta, p_2 + mu, omegap, &g2g);
	//stencil_display(&g2g, 1.0);

	stencil_naive(p, a, h, p_2+mu, omegap, k, c, &g2g);

	//	Pre-processing (Top level)
	stencil_initialize(&gamma, (long) ceil(2.0*alpha)+p_2+mu, STENCIL_SHAPE_CUBE);
	stencil_populate(&gamma, c, k, STENCIL_FUNCTION_TYPE_GAMMA, h/a);
	//stencil_display(&gamma, h/a);

	stencil_initialize(&tg2g, gamma.size-p_2-mu, gamma.shape);
	stencil_shift_infinite(&gamma, p_2+mu, omegap, &tg2g);
	//stencil_display(&tg2g, 1.0);

	stencil_naive_top(p, a, h, p_2+mu, omegap, k, c, &tg2g);

	//	Free dynamically allocated memory
	dynfree(c);
	dynfree(omegap);
	stencil_free(&theta);
	stencil_free(&gamma);
	stencil_free(&g2g);
	stencil_free(&tg2g);
*/
}

void test_softener(SOFTENER* s)
{
	long		Samples = 100;
	long		Len = Samples + 1;
	long		i = 0;
	double*		X = NULL;
	double*		FX = NULL;
	double*		DFX = NULL;
	double		fx = 0.0;
	double		dfx = 0.0;
	void		(*MyGamma)(void*,long,double*,double*,double*);
	void		(*MyTheta)(void*,long,double*,double*,double*);
	void*		Init = NULL;
	SOFTENER*	t = NULL;

	if (s == NULL)
	{
		Init = &even_powers_initialize;
		softener_initialize((void**)&t, sizeof(EVEN_POWERS), Init, 4);
		s = t;
	}

	X = (double*) dynvec(Len, sizeof(double));
	FX = (double*) dynvec(Len, sizeof(double));
	DFX = (double*) dynvec(Len, sizeof(double));

	for (i = 0; i < Len; i++)
	{
		X[i] = (double) 2.5*i / (double) Samples;
	}

	MyTheta = s->split;
	(*MyTheta)(s, Len, X, FX, DFX);

	for (i = 0; i < Len; i++)
	{
		fx = theta(s->p2p, s->k, X[i], &dfx);
		printf("i=%03ld, x=%f |f-F| = %e, |df-DF| = %e\n", i, X[i], fabs(fx-FX[i]), fabs(dfx-DFX[i]));
	}

	printf("\n");
	MyGamma = s->soften;
	(*MyGamma)(s, Len, X, FX, DFX);

	for (i = 0; i < Len; i++)
	{
		fx = _gamma(s->p2p, s->k, X[i], &dfx);
		printf("i=%03ld, x=%f |f-F| = %e, |df-DF| = %e\n", i, X[i], fabs(fx-FX[i]), fabs(dfx-DFX[i]));
	}

	if (t != NULL)
	{
		t->uninitialize(t);
	}

	dynfree(X);
	dynfree(FX);
	dynfree(DFX);
}
#endif

//	End of file