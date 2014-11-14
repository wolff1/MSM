//-------|---------|---------|---------|---------|---------|---------|---------|
/*
tester.h - 
*/

#ifndef TESTER_H
#define TESTER_H

#include "all.h"
#include "b_spline.h"
#include "c1_spline.h"
#include "even_powers.h"
#include "interpolant.h"
#include "memory.h"
#include "msm.h"
#include "method.h"
#include "output.h"
#include "polynomial.h"
#include "softener.h"
#include "stencil.h"

//	**** POLYNOMIAL ****
void test_mpoly(void);

//	**** INTERPOLANT ****
void phi_test_all(void);
void print_nesting_coefficients(void);
void phi_nesting_test(void);

//	**** B_SPLINE ****
void test_blurring_operator(void);
void test_omega(void);
void test_omegap(void);
void test_convert_to_shifts(void);

//	**** C1_SPLINE ****
void driverC1(void);
void phi_test_allC1(void);
void print_nesting_coefficientsC1(void);
void phi_nesting_testC1(void);

//	**** EVEN_POWERS ****
void test_thetas(void);
void splitting_test(void);
void plot_splittings(int samples, int nlev, double a, double d,
						double* X, double** F);
void gamma_test_all(void);

//	**** OTHER ****
void test_sinc(void);
void test_preprocessing(void);

#endif

//	End of file