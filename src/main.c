//-------|---------|---------|---------|---------|---------|---------|---------|
/*
main.c - entry point for the tester
*/

#include "tester.h"

void test_mkl_MMM(void);
void mkl_memory_check(void);
void test_sinc(void);

void msm_init(void);
void msm_preprocess(void);
void msm_eval(void);

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
		printf("*15 - Test Preprocessing         *\n");
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

			case 15:
				test_preprocessing();
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
	compute_omega(4, omega_max, omega3);
	compute_omega(6, omega_max, omega5);
	compute_omega(8, omega_max, omega7);

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

void msm_init(void)
{
	//	input: accuracy parameters:	a, h, p, mu, k, stencil shape,
	//			other parameters:	L (domain size), lmax (number of grids)
	//			internal objects:	g2p, g2fg, p2p, omega', g2g, tg2g
	//			grids:				grid, 1 to lmax

	//	g2p (g2p is p/2 x p)
	//compute_g2p(short p, double* g2p);

	//	g2fg (len = p/2+1)
	//compute_g2fg(short p, double* g2fg);

	//	p2p (p2p is k+1)
	//gamma_init(short k, double* p2p);

	//	omega' (omegap is p/2+mu+1)
	//compute_omega_prime(short p, short mu, double* omegap);

	//	g2g
	//	tg2g
	msm_preprocess();
}

void msm_preprocess(void)
{
	//	input: omegap', alpha = a/h, L = domain size, K-top, p, p_2, mu

	//	if g2g == NULL
	//		compute intermediate stencil (only happens once)
	//stencil_initialize(&theta, (long) ceil(2.0*alpha), STENCIL_SHAPE_SPHERE);
	//stencil_populate(&theta, c, k, STENCIL_FUNCTION_TYPE_THETA, h/a);
	//stencil_initialize(&g2g, theta.size, theta.shape);
	//stencil_shift(&theta, p_2 + mu, omegap, &g2g);
	//	write theta to a file for future use and to free up some memory?

	//	if tg2g == NULL or domain enlarged
	//		first full computation is resizing from 0 to X
	//		additional computations are resizing from X to Y
	//stencil_initialize(&gamma, 10 * (long) ceil(2.0*alpha), STENCIL_SHAPE_CUBE);
	//stencil_populate(&gamma, c, k, STENCIL_FUNCTION_TYPE_GAMMA, h/a);
	//stencil_initialize(&tg2g, gamma.size, gamma.shape);
	//stencil_shift(&gamma, p_2 + mu, omegap, &tg2g);
	//	write gamma to a file for future use and to free up some memory?
}

void msm_eval(void)
{
	//	input: r, q, options, parameters

	//	short_range
	//	anterpolate
	//		direct
	//		restrict
	//	direct_top
	//		prolongate
	//	interpolate
	//	exclusions
}

// End of file
