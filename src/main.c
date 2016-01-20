//-------|---------|---------|---------|---------|---------|---------|---------|
/*
main.c - entry point for the tester
*/
#include <time.h>

#include "msm.h"
#include "simulator.h"
#include "tester.h"

void test_mkl_MMM(void);

void test_msm_preprocessing(void);
void test_simulator(void);
void test_parallel_division(void);
void test_simulator_water(void);
void single_splitting(void);
double weighted_fn(void);
void test_quasi_interp_1d(void);

/*
#include "tester.h"
int main(int argc, char* argv[])
{
	gamma_test_all();
	return 0;
}
*/

int main(int argc, char* argv[])
{
	int				choice = 1;

	while (choice > 0)
	{
		printf("************** MENU **************\n");
		printf("* 1 - Test MSM Preprocessing     *\n");
		printf("* 2 - Simulator                  *\n");
		printf("* 3 - Parallel division of procs *\n");
		printf("* 4 - Simulation Water(s)        *\n");
		printf("*                                *\n");
		printf("* 5 - Figure 1 (Splitting)       *\n");
		printf("* 6 - Figure 2 (Sinc)            *\n");
		printf("* 7 - Figure 4 (B-Spline Nesting)*\n");
		printf("* 8 - Figure 5 (C1 Not Nesting)  *\n");
		printf("*                                *\n");
		printf("* 9 - Single splitting (1D)      *\n");
		printf("*10 - Quasi-Interp (1D)          *\n");
		printf("**********************************\n");
		printf("* 0 - Exit                       *\n");
		printf("**********************************\n");
		printf("Your selection: ");
		scanf("%d", &choice);

		// Give them what they want.
		switch (choice)
		{
			case 1:	//	Test the preprocessing routines
				test_msm_preprocessing();
				break;

			case 2:	//	Test the simulator
				test_simulator();
				break;

			case 3:
				test_parallel_division();
				break;

			case 4:
				test_simulator_water();
				break;

			case 5:	// Figure 1
				splitting_test();
				break;

			case 6:	//	Figure 2
				test_sinc();
				break;

			case 7:	//	Figure 4
				phi_nesting_test();
				break;

			case 8:	//	Figure 5
				phi_nesting_testC1();
				break;

			case 9:	//	Single-splitting (1D)
				single_splitting();
				break;

			case 10://	quasi-interpolation (1d)
				test_quasi_interp_1d();
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
int main(int argc, char* argv[])
{
	int				choice = 1;

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
		printf("*16 - Test Simulator             *\n");
		printf("*17 - Test Softening (vector)    *\n");
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

			case 16:
				simulator();
				break;

			case 17:
				test_softener(NULL);

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
*/

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

void test_msm_preprocessing(void)
{
	//	Generic Variables
	void*		Ptr = NULL;
	METHOD*		Method = NULL;
	size_t		Size = 0;
	void*		Init = NULL;

	if (1)
	{	//	MSM
		Size = sizeof(MSM);
		Init = &msm_initialize;
	}

	//	Initialize method
	Ptr = (METHOD*) dynmem(Size);
	method_initialize(Ptr, Init, 0, NULL, NULL);	//	FIXME: 2 NULLs will not work!
	Method = (METHOD*) Ptr;

	//	Run simulation
	(*Method->preprocess)(Method, 10.0);
	//(*Method->evaluate)(Method, 0, NULL, NULL, NULL);
	(*Method->uninitialize)(Method);
	dynfree(Method);
}

void test_simulator(void)
{
	SIMULATOR*	MySimulator = (SIMULATOR*) dynmem(sizeof(SIMULATOR));

	simulator_initialize(MySimulator);
	simulator_run(MySimulator);
	simulator_examine_results(MySimulator);
	simulator_uninitialize(MySimulator);

	dynfree(MySimulator);
}

void test_parallel_division(void)
{
	int		xp = 0;
	int		yp = 0;
	int		zp = 0;
	int		x = 0;
	int		y = 0;
	int		z = 0;
	int		q = 0;
	int		m = 0;
	int		n = 0;
	int		idx = 0;

	for (zp = -2; zp <= 2; zp++)
	{
		for (yp = -2; yp <= 2; yp++)
		{
			for (xp = -2; xp <= 2; xp++)
			{
//				q = MIN(abs(x),MIN(abs(y),abs(z)));
//				n = MAX(abs(x),MAX(abs(y),abs(z)));
//				m = MIN(abs(x)-q,abs(y)-q) + MIN(abs(y)-q,abs(z)-q) + MIN(abs(x)-q,abs(z)-q) + q;

//				z = abs(zp);	y = abs(yp);	x = abs(xp);
//				if (x < y)	{if (x < z) {q = x;	if (y < z) {m = y; n = z;} else {m = z; n = y;}} else {q = z; m = x; n = y;}}
//				else		{if (y < z) {q = y;	if (x < z) {m = x; n = z;} else {m = z; n = x;}} else {q = z; m = y; n = x;}}
//				idx = STENCIL_MAP_X(q) + STENCIL_MAP_Y(m) + STENCIL_MAP_Z(n);

				SORT_RTN_STNCL_IDX(abs(xp),abs(yp),abs(zp),q,m,n,idx);
				printf("(%+d,%+d,%+d) -> (%+d,%+d,%+d) -> (%+d,%+d,%+d) -> %+d\n", xp,yp,zp, abs(x),abs(y),abs(z), q,m,n, idx);
			}
		}
	}
/*
	long		N = 0;				//	# particles
	long		K = 0;				//	# grid points
	long		P = 0;				//	# procs
	double		a = 0.0;			//	cut-off distance
	double		h = 0.0;			//	grid spacing
	short		p = 0;				//	Order of interpolant

	double		alpha = 0.0;		//	Half stencil radius
	double		Ps = 0.0;			//	# procs split 1
	double		Pt = 0.0;			//	# procs split 2
	long		D0 = 0;
	long		D1 = 0;
	long		D2 = 0;
	long		R1 = 0;
	long		R2 = 0;
	long		P1 = 0;
	long		P0 = 0;
	double		Kfake = 0.0;
	double		Kconst = 0.0;

	printf("Please enter particles N: ");
	scanf("%ld", &N);
	printf("Please enter grid points K: ");
	scanf("%ld", &K);
	printf("Please enter processors P: ");
	scanf("%ld", &P);
	printf("Please enter cut-off a: ");
	scanf("%lf", &a);
	printf("Please enter grid spacing h: ");
	scanf("%lf", &h);
	printf("Please enter order p: ");
	scanf("%hd", &p);

	alpha = 2.0*a/h;
	printf("2*alpha = %f\n", alpha);
*/
/*
	Work per step = (Work per point)*(number of points)
	Time per step = (Time/Work)*Work
		=> 1 processor can do 1 unit Work in 1 unit time
		=> We're not trying to balance work, per se, we're trying to balance *time*

	D0 ~= R1 + D1 + P0:
	[(32*pi/3)*(alpha)^3]*K ~= [(p+1)^3]*K + [(4*pi/3)*(alpha)^3]*K + [(p+1)^3]*K
	[(32*pi/3)*(alpha)^3]    = [(p+1)^3] +   [(4*pi/3)*(alpha)^3] +   [(p+1)^3]

	T(D0,P-S) = T(R1,S) + T(D1,S-T) + T(P0, S):
	[(32*pi/3)*(alpha)^3]/(P-S)  =   [(p+1)^3]/S +   [(4*pi/3)*(alpha)^3]/(S-T) +   [(p+1)^3]/S
	[(32*pi/3)*(alpha)^3]/(P-S)  = 2*[(p+1)^3]/S +   [(4*pi/3)*(alpha)^3]/(S-T)

	//	D{0} to D{L-2}
	T(D{i},P{i}) = T(R{i+1},P{i-1}-P{i}) + T(D{i+1},P{i}-P{i+1}) + T(Pr{0}, P{i}):
	[(32*pi/3)*(alpha)^3]/(P{i})  = 2*[(p+1)^3]/(P{i-1}-P{i}) +   [(4*pi/3)*(alpha)^3]/(P{i}-P{i+1}); 0 <= i <= L-3; P{-1} = P; 

	D{i}   = (32*pi/3)*(alpha)^3;
	R{i+1} = (p+1)^3;
	D{i+1} = (4*pi/3)*(alpha)^3;	//	-> (1/8)*D{i}
	Pr{i}  = (p+1)^3;
	D{top} = 8^{-(L-2)}K

	P{i}/((32*pi/3)*(alpha)^3)  = (P{i-1}-P{i})/(2*(p+1)^3) +   (P{i}-P{i+1})/((4*pi/3)*(alpha)^3)
	P{i}/(D{i})  = (P{i-1}-P{i})/(R{i+1}+Pr{i}) +   (P{i}-P{i+1})/(D{i+1})
*/
}

void test_simulator_water(void)
{
	SIMULATOR*	MySimulator = (SIMULATOR*) dynmem(sizeof(SIMULATOR));

	simulator_initialize(MySimulator);
	simulator_run_water(MySimulator);
//	simulator_examine_results(MySimulator);
	simulator_uninitialize(MySimulator);

	dynfree(MySimulator);
}

void single_splitting(void)
{
	long			i = 0;
	long			j = 0;
	long			samples = 10;
	long			M = 0;
	double			err = 0.0;
	double			num = 1.0;
	double			den = 1.0;
	double			x1 = 0.0;
	double			x2 = 1.0;
	double			f1 = 0.0;
	double			f2 = 0.0;
	double*			X = NULL;
	double*			F = NULL;
	FILE*			fp = NULL;
	char			filename[256];
	time_t			rawtime;
	struct tm*		now;
	double			min_D = 0.0;

	MSM_PARAMETERS	mp;
	B_SPLINE*		bs = NULL;
	EVEN_POWERS*	ep = NULL;

	//	Set up necessary parameters
	//mp.h = 2.5;
	//mp.a = 4.0*mp.h;
	//mp.mu = 10;
	//mp.p = 4;
	//mp.k = mp.p;	// k = p/2 - 1 should give EXACT
	//mp.L = 4;
	//mp.D = 10;
	msm_parameters_input(&mp);
	printf("Samples (Default: %ld): ", samples);
	scanf("%ld", &samples);

	//	Create B-Spline object
	bs = (B_SPLINE*) dynmem(sizeof(B_SPLINE));
	interpolant_initialize(bs, (void*)b_spline_initialize, &mp);

	//	Create Gamma object
	ep = (EVEN_POWERS*) dynmem(sizeof(EVEN_POWERS));
	softener_initialize(ep, (void*)even_powers_initialize, mp.k);

	//	Output:		expected error based on method parameters & analysis (from slide 6/13 of 10/01/2015 thesis presentation)
	for (i = 1; i <= mp.p/2; i++)
	{
		err += 1.0/(double)i;
		err += 2.0/(double)(i+mp.p/2.0);
	}
	err -= log(2.0*mp.p);
//printf("gamma_%d = %f\n", mp.p, err);
	err += 2.0*log(4.0/PI);
	err += log(mp.p-1);
	err *= 2.0/PI;
	for (i = 2; i <= mp.p; i+=2)
	{
		num *= (double)i-1.0;
		den *= (double)i;
	}
//printf("factorial thing: %f\n", num/den);
	err += (num/den);
	err *= pow(mp.h/2.0, mp.p);
	ep->cmn.derivative(ep,1,&x1,&f1,mp.p);
	ep->cmn.derivative(ep,1,&x2,&f2,mp.p);
//printf("gamma^{(%d)}(%f) = %f, gamma^{(%d)}(%f) = %f\n", mp.p, x1, f1, mp.p, x2, f2);
	err *= MAX(abs(f1),abs(f2));
	err /= pow(mp.a, mp.p+1);
printf("expected error is %e\n", err);

	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));

	//	Get current time and format 
	time(&rawtime);
	now = localtime(&rawtime);
/*
	//	Open file for writing and then write stuff
	strftime(filename, sizeof(filename), "../../../data/dgamma_%Y-%m-%d_%H-%M-%S.dat", now);
	if ((fp = fopen(filename, "w")) != NULL)
	{
		fprintf(fp, "# This is: %s, derivatives of gamma\n", filename);
		fprintf(fp, "# Parameters: \n");
		fprintf(fp, "#\ta = %lf\n", mp.a);
		fprintf(fp, "#\th = %lf\n", mp.h);
		fprintf(fp, "#\tp = %hd\n", mp.p);
		fprintf(fp, "#\tk = %hd\n", mp.k);
		fprintf(fp, "#\tmu = %hd\n", mp.mu);
		fprintf(fp, "#\tD = %lf\n", mp.D);
		fprintf(fp, "#\tL = %hd\n\n", mp.L);

		//	Output "high" derivative of gamma over domain
		for (i = 0; i <= samples; i++)
		{
			X[i] = 1.0 + ((double)i/(double)samples)*mp.D;
		}

		for (i = 0; i <= mp.k; i++)
		{
			fprintf(fp, "# %ld-th derivative\n", i);

			//	Compute i-th derivatie for all "samples" over domain
			ep->cmn.derivative(ep, samples+1, X, F, i);

			//	write derivative(s) to file
			for (j = 0; j <= samples; j++)
			{
				fprintf(fp, "%04ld\t%e\t%e\n", j, X[j], F[j]);
			}
			fprintf(fp, "\n");
		}

		//	Close file
		fclose(fp);
		printf("*** COMPLETED GAMMA DERIVATIVES TEST. Results in file: %s ***\n", filename);
	}
*/
	min_D = mp.h*pow(2,mp.L-1)*1.0;
	//	Set up number of fine-grid points
	if (min_D > mp.D)
	{
		printf("*** Updating D from %lf to ", mp.D);
		mp.D = min_D;
		printf("%lf to ensure coarse grid is big enough. ***\n", mp.D);
	}
	M = floor(mp.D / mp.h) + 1;	// # of grid points necessary for domain (but NOT "expansion")

	for (i = mp.L; i > 0; i--)
	{
		printf("|Grid %ld| = %ld\n", i, (M-1)/(long)pow(2.0, i-1) + mp.p + 1);
	}

//	display_vector_d(bs->omega, bs->cmn.p/2 + bs->mu + 1);
	dynfree(bs->omega);
	b_spline_compute_omega(bs, 100);
//	display_vector_d(bs->omega, 100);

	//	Level L:	compute gamma
	//				create interpolant of gamma on grid L
	//				compute values of interpolant over domain

	//	Level L-1:	compute gamma - interp(L)
	//				create interpolant of gamma - interp(L) on grid L-1
	//				compute values of interpolant over domain

	//	Level ell:	compute gamma - interp(ell+1)
	//				create interpolant of gamma - interp(ell+1) on grid ell
	//				compute values of interpolant over domain

	//	Level 1:	compute gamma - interp(2)
	//				create interpolant of gamma - interp(2) on grid 1
	//				compute values of interpolant over domain

	//	Output:		gamma, interp(L) + ... + interp(1), interp(L), interp(L-1), ..., interp(2), interp(1)

	//	Destory X and F vectors
	dynfree(X);
	dynfree(F);

	//	Destroy B-Spline object
	bs->cmn.uninitialize(bs);
	dynfree(bs);

	//	Destroy Gamma object
	ep->cmn.uninitialize(ep);
	dynfree(ep);
}

double weighted_fn(void)
{
	return 0.0;
}

void b_spline_compute_omega_prime_1d(short p, short mu, double* omegap)
{
	long			i = 0;
	double*			B = NULL;
    double*			BE = NULL;
	double*			Ap = NULL; // A' ~= B^{-1}
	double*			Ctmp = NULL;
	double*			C = NULL;
	double*			CE = NULL;
    double*			c = NULL;
	double*			pd = NULL;
	double*			pE = NULL;
	double*			c_full = NULL;
	double*			p_full = NULL;
	double*			high_terms = NULL;
	short			p_2 = p/2;

	//	Build B
	B = (double*) dynvec(p_2, sizeof(double));
	b_spline_compute_blurring_operator(p_2-1, B);

	//	Build A'
	Ap = (double*) dynvec(p_2, sizeof(double));
	b_spline_compute_operator_inverse(p_2-1, B, p_2-1, Ap);

	//	Compute C
	Ctmp = (double*) dynvec(p-1, sizeof(double));
	mpoly(p_2-1, B, p_2-1, Ap, Ctmp);
	//	keep only (delta^2)^{p/2} and higher
	C = (double*) dynvec(p_2-1, sizeof(double));
	for (i = p_2; i < p; i++)
	{
		//C[i-p_2] = -Ctmp[i]; // -Ctmp?
		C[i-p_2] = Ctmp[i];
	}
	//	convert delta^2 to shifts
	CE = (double*) dynvec(p_2-1, sizeof(double));
	b_spline_convert_delta2_to_shifts(p_2-2, C, CE);
	b_spline_convert_delta2_to_shifts(p_2-1, Ap, omegap);

	//	Solve for E in shift operators
    BE = (double*) dynvec(p_2, sizeof(double));
    b_spline_convert_delta2_to_shifts(p_2-1, B, BE);
    c = (double*) dynvec(mu+1, sizeof(double));
    bibst_lss(-1, pow(2.0,-53), p_2, BE, p_2-1, p_2-1, CE, mu+1, c);

	//	Convert (delta^2)^{p/2} to shifts
	pd = (double*) dynvec(p_2+1,sizeof(double));
	pE = (double*) dynvec(p_2+1,sizeof(double));
	pd[p_2] = 1.0;
	b_spline_convert_delta2_to_shifts(p_2,pd,pE);

	//	Form symmetric pE and c operators
	c_full = (double*) dynvec(2*mu+1,sizeof(double));
	c_full[mu] = c[0];
	for (i = 1; i <= mu; i++)
	{
		c_full[mu+i] = c[i];
		c_full[mu-i] = c[i];
	}

	p_full = (double*) dynvec(p+1,sizeof(double));
	p_full[p_2] = pE[0];
	for (i = 1; i <= p_2; i++)
	{
		p_full[p_2+i] = pE[i];
		p_full[p_2-i] = pE[i];
	}

	//	Apply operators to one another
	high_terms = (double*) dynvec(p+2*mu+1,sizeof(double));
	mpoly(2*mu, c_full, p, p_full, high_terms);

	//	Compute A = A'(E) + E(E)
	for (i = 0; i <= p_2+mu; i++)
	{
		omegap[i] += high_terms[p_2+mu+i];
	}

	// Free dynamically allocated memory
	dynfree(B);
    dynfree(BE);
	dynfree(Ap);
	dynfree(Ctmp);
	dynfree(C);
	dynfree(CE);
    dynfree(c);
	dynfree(pd);
	dynfree(pE);
	dynfree(c_full);
	dynfree(p_full);
	dynfree(high_terms);
}

void setup_interpolation_grid(short p, short mu, long nodes, double* G, double* Gmin, double* Gmax)
{
	long	i = 0;
	double	t1 = 0.0;
	double	t2 = 0.0;

	//	This could be a prompt for user input
	if (((*Gmin) == 0.0) && ((*Gmax) == 0.0))
	{
		t1 = (double) rand();
		t2 = (double) rand();
		(*Gmin) = MIN(t1, t2);
		(*Gmax) = MAX(t1, t2);
	}

	printf("Domain: [%lf, %lf]\n\n", (*Gmin), (*Gmax));

	//	Size of G is nodes+1+2*(p+mu-1)
	//	nodes+1 between Gmin and Gmax
	//	The rest are before Gmin and after Gmax

	for (i = -(p+mu-1); i <= nodes+p+mu-1; i++)
	{
		G[p+mu-1+i] = (*Gmin) + ((*Gmax)-(*Gmin))*(double)i/(double)nodes;
	}
}

double eval_poly(short p, double* poly, double x)
{
	double	fx = 0.0;
	long	j = 0;

	//	Horner's Rule:
	fx = poly[p-1];
	for (j = p-2; j >= 0; j--)
	{
		fx = fx*x + poly[j];
	}

	return fx;
}

void compute_ordinates_on_interpolation_grid(short p, long Len, double* G, double* F, double* poly, EVEN_POWERS* ep)
{
	long			i = 0;
	long			j = 0;
	double*			DF = NULL;

	if (ep == NULL)
	{
		printf("Polynomial is: \np(x) = ");
		for (i = 0; i < p; i++)
		{
			poly[i] = (double) rand() / (double) rand();
			printf("%lf*x^%ld %s ", poly[i], i, (i < p-1 ? "+" : "\n"));
		}
		printf("\n");

		for (i = 0; i < Len; i++)
		{
			F[i] = eval_poly(p, poly, G[i]);
		}
	}
	else
	{
		DF = (double*) dynvec(Len, sizeof(double));
		ep->cmn.soften(ep, Len, G, F, DF);
		dynfree(DF);
	}
}

void compute_interpolant(short p, short mu, long nodes, double* omegap, double* F, double* I)
{
	long			i = 0;
	long			j = 0;
	short			p_2 = p/2;

	for (j = 0; j <= nodes+2*(p_2-1); j++)
	{
		I[j] += omegap[i]*F[p_2+mu+j];
	}

	for (i = 1; i <= p_2+mu; i++)
	{
		//	i represents shift: E^i
		for (j = 0; j <= nodes+2*(p_2-1); j++)
		{
			I[j] += omegap[i]*(F[p_2+mu+j+i] + F[p_2+mu+j-i]);
		}
	}
}

void interpolate_over_domain(short p, double h, double Gmin, double* I, long samples, double* X_BAR, double* F_BAR)
{
	long			i = 0;
	long			j = 0;
	short			p_2 = p/2;
	MSM_PARAMETERS	mp;
	B_SPLINE*		bs = NULL;
	double			xmin = 0.0;
	double*			BSX = NULL;
	double*			BSF = NULL;
	double*			DBSF = NULL;

//	msm_parameters_input(&mp);
	mp.a = 1.0;
	mp.alpha = 1.0;
	mp.h = h;
	mp.D = X_BAR[samples] - X_BAR[0];
	mp.L = 1;
	mp.mu = 10;
	mp.p = p;
	mp.k = mp.p;

	//	Create B-Spline object
	bs = (B_SPLINE*) dynmem(sizeof(B_SPLINE));
	interpolant_initialize(bs, (void*)b_spline_initialize, &mp);

	BSX = (double*) dynvec(p, sizeof(double));
	BSF = (double*) dynvec(p, sizeof(double));
	DBSF = (double*) dynvec(p, sizeof(double));

	//	f_bar(x) = \sum f_hat(m) * Phi (x/h - m)
	for (i = 0; i <= samples; i++)
	{
		//	Set up input for BS evaluation
		xmin = floor((X_BAR[i]-Gmin)/mp.h) - p_2 + 1;
		for (j = 0; j < p; j++)
		{
			//	x/h-m
			BSX[j] = (X_BAR[i]-Gmin)/mp.h - (xmin + (double) j);
		}

		//	Evaluate BS at BSX
		b_spline_evaluate(bs, p, BSX, BSF, DBSF);

		//	Use BSF and I to interpolate value for G[i]
		F_BAR[i] = 0.0;
		for (j = 0; j < p; j++)
		{
			F_BAR[i] += I[p_2-1+(long)xmin+j]*BSF[j];
		}
	}

	//	Destroy B-Spline object
	bs->cmn.uninitialize(bs);
	dynfree(bs);

	dynfree(BSX);
	dynfree(BSF);
	dynfree(DBSF);
}

void test_quasi_interp_1d(void)
{
	short			p = 0;
	short			p_2 = 0;
	long			i = 0;
	long			j = 0;
	short			mu = 1;

	double*			omegap = NULL;
	double			Gmin = 0.0;
	double			Gmax = 0.0;

	time_t			t;
	long			nodes = 10;
	long			samples = 10;
	double*			poly = NULL;
	double*			G = NULL;
	double*			F = NULL;
	double*			I = NULL;

	EVEN_POWERS*	ep = NULL;

	double*			X_BAR = NULL;
	double*			F_BAR = NULL;
	double*			F_HAT = NULL;
	double*			DF_HAT = NULL;
	double*			DIFF = NULL;
	double*			RDIFF = NULL;

	double			l1 = 0.0;
	double			l2 = 0.0;
	double			loo = 0.0;

	FILE*			fp = NULL;
	char			fn[64];

//	GET NECESSARY PARAMETERS
	srand((unsigned) time(&t));
	printf("p = ");
	scanf("%hd", &p);
	p_2 = p/2;

	printf("mu = ");
	scanf("%hd", &mu);

	printf("Interpolation Nodes (Default: %ld): ", nodes);
	scanf("%ld", &nodes);

//	BUILD INTERPOLATION OPERATOR 1D
	omegap = (double*) dynvec(p_2+mu+1, sizeof(double));
	b_spline_compute_omega_prime_1d(p, mu, omegap);

//	COMPUTE DISCRETE FUNCTION VALUES
	//	Create sequence of values, fx, for x in [min, max]
	G = (double*) dynvec(nodes+1+2*(p+mu-1), sizeof(double));
	printf("Gmin = ");
	scanf("%lf", &Gmin);
	printf("Gmax = ");
	scanf("%lf", &Gmax);
	setup_interpolation_grid(p, mu, nodes, G, &Gmin, &Gmax);

	//	Build polynomial of degree p-1, for which our interpolant *should* be exact
	F = (double*) dynvec(nodes+1+2*(p+mu-1), sizeof(double));

	//	Create Gamma object and/or polynomial object
	ep = (EVEN_POWERS*) dynmem(sizeof(EVEN_POWERS));
	softener_initialize(ep, (void*)even_powers_initialize, p/2-1);	// polynomial part of spline with degree < p-1 yields better accuracy
	poly = (double*) dynvec(p, sizeof(double));

	compute_ordinates_on_interpolation_grid(p, nodes+1+2*(p+mu-1), G, F, poly, ep);

//	APPLY INTERPOLATION OPERATOR (i.e. anti-blur function values)
	I = (double*) dynvec(nodes+1+2*(p_2-1), sizeof(double));
	compute_interpolant(p, mu, nodes, omegap, F, I);

//	INTERPOLATE AND CHECK ACCURACY
	printf("# of Samples (Default: %ld): ", samples);
	scanf("%ld", &samples);

	X_BAR = (double*) dynvec(samples+1, sizeof(double));
	F_BAR = (double*) dynvec(samples+1, sizeof(double));
	F_HAT = (double*) dynvec(samples+1, sizeof(double));
	DIFF = (double*) dynvec(samples+1, sizeof(double));
	RDIFF = (double*) dynvec(samples+1, sizeof(double));

	//	Set up X_BAR
	for (i = 0; i <= samples; i++)
	{
		X_BAR[i] = Gmin + (Gmax-Gmin)*i/samples;
	}

	interpolate_over_domain(p, G[1]-G[0], Gmin, I, samples, X_BAR, F_BAR);

	DF_HAT = (double*) dynvec(samples+1, sizeof(double));
	if (ep == NULL)
	{
		//	Compute target function
		for (i = 0; i <= samples; i++)
		{
			F_HAT[i] = eval_poly(p, poly, X_BAR[i]);
		}
	}
	else
	{
		ep->cmn.soften(ep, samples+1, X_BAR, F_HAT, DF_HAT);
	}

	//	Compute error(s)
	for (i = 0; i <= samples; i++)
	{
		DIFF[i] = fabs(F_HAT[i] - F_BAR[i]);
		RDIFF[i] = DIFF[i] / fabs(F_HAT[i]);
	}

	//	Display norm(s) of error
	j = 0;
	for (i = 0; i <= samples; i++)
	{
		j++;
		l1 += fabs(RDIFF[i]);
		l2 += fabs(RDIFF[i])*fabs(RDIFF[i]);
		loo = MAX(loo, RDIFF[i]);
	}
	l2 = sqrt(l2);
	printf("j=%ld, L1=%e, L2=%e, Loo=%e\n", j, l1, l2, loo);

	//	Write output file to plot via gnuplot
	strftime(fn, 64*sizeof(char), "1DInterp_%Y_%m_%d_%H_%M_%S.dat", localtime(&t));
	if ((fp = fopen(fn, "w")) != NULL)
	{
		//	Write samples to file
		for (i = 0; i <= samples; i++)
		{
			//	Output: i, X_BAR, F_HAT, F_BAR, DIFF, RDIFF
			fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, X_BAR[i], F_HAT[i], F_BAR[i], DIFF[i], RDIFF[i]);
		}

		//	Close file
		fclose(fp);
		printf("Output file written!\n");
	}

	//	Free dynamically allocated memory
	dynfree(poly);
	dynfree(G);
	dynfree(F);
	dynfree(I);

	dynfree(X_BAR);
	dynfree(F_BAR);
	dynfree(F_HAT);
	dynfree(DF_HAT);
	dynfree(DIFF);
	dynfree(RDIFF);

	//	Destroy Gamma object
	ep->cmn.uninitialize(ep);
	dynfree(ep);

	dynfree(omegap);	//	REMOVE THIS FOR REAL CODE!
}

//void test_parallel_division(void)
//{
//	long		N = 0;				//	# particles
//	long		K = 0;				//	# grid points
//	long		P = 0;				//	# procs
//	double		Ps = 0.0;			//	# procs split 1
//	double		Pt = 0.0;			//	# procs split 2
//	double		a = 0.0;			//	cut-off distance
//	double		h = 0.0;			//	grid spacing
//	double		alpha = 0.0;		//	Half stencil radius
//	short		p = 0;				//	Order of interpolant
//	long		D0 = 0;
//	long		D1 = 0;
//	long		D2 = 0;
//	long		R1 = 0;
//	long		R2 = 0;
//	long		P1 = 0;
//	long		P0 = 0;
//	double		Kfake = 0.0;
//	double		Kconst = 0.0;
//
//	printf("Please enter particles N: ");
//	scanf("%ld", &N);
////	printf("Please enter grid points K: ");
////	scanf("%ld", &K);
//	K = N;
//	printf("Please enter processors P: ");
//	scanf("%ld", &P);
//	P = MAX(P,100);
//	printf("Please enter cut-off a: ");
//	scanf("%lf", &a);
//	printf("Please enter grid spacing h: ");
//	scanf("%lf", &h);
//	alpha = 2.0*a/h;
//	printf("Please enter order p: ");
//	scanf("%hd", &p);
//	printf("2*alpha = %f\n", alpha);
//
////	printf("You entered: N=%ld, K=%ld, P=%ld\n", N, K, P);
//
//	//	4-level method
//	D0 = (long)ceil(pow(alpha,3))*K;
//	D1 = (long)ceil(pow(alpha,3))*(long)ceil((double)K/8.0);
//	D2 = (long)ceil((double)K*K/pow(64.0,2));
//	R1 = p*p*p*(long)ceil(K/pow(8.0,1));
//	R2 = p*p*p*(long)ceil(K/pow(8.0,2));
//	P1 = R2;
//	P0 = R1;
//	printf("Here are the complexities:\nD0=%ld\nR1=%ld, D1=%ld, P0=%ld\nR2=%ld, D2=%ld, P1=%ld\n", D0, R1, D1, P0, R2, D2, P1);
//
//	srand(time(NULL));
//	Kconst = (double)rand()/((double)RAND_MAX);
//	Kconst = 0.05 + (0.20)*Kconst;
//	Kfake = (double)D0 + Kconst*D0;
////	printf("Kfake = %f, Kconst=%f\n", Kfake, Kconst);
////	printf("D0/Kfake = %f\n", (double)D0/Kfake);
//	Ps = (double)P - (double)P*(double)D0/Kfake;		//	D0/K should be close to 1 so that Ps is relatively small
//	Pt = Ps*(R2+D2+P1)/(double)(D1+R2+D2+P1);	//	Pt will be even less than Ps
////	printf("P = %ld, s = %f, t = %f, Ps = %ld, Pt = %ld\n", P, Ps, Pt, (long)ceil(Ps), (long)ceil(Pt));
////	printf("\n\n");
//	Ps = ceil(Ps);
//	Pt = ceil(Pt);
//	printf("P = %ld, P-s = %ld, s = %ld, s-t = %ld, t = %ld\n", P, P-(long)Ps, (long)Ps, (long)Ps-(long)Pt, (long)Pt);
//
//	printf("#1: D0 / (P-s) = K / P, %f ?= %f\n", (double)D0/((double)P-Ps), Kfake / (double)P);
//	printf("#2: D0 / (P-s) = R1/s + D1/(s-t) + P0/s => %f ?= %f\n", (double)D0/((double)P-Ps), (double)R1/Ps + (double)D1/(Ps-Pt) + (double)P0/Ps);
//	printf("#3: D1 / (s-t) = R2/t + D2/t + P1/t => %f ?= %f\n", (double)D1/(Ps-Pt), ((double)(R2 + D2 + P1))/Pt);
//}

// End of file