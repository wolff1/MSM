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

			case 5:	// This will produce plot (use OSX) for smoothing
				splitting_test();
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