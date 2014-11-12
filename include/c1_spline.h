//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.h - C1 phi
*/

#define PHI_DATA_C1		"phiC1.dat"

#define PHI_DATA_C1_LEN	strlen(PHI_DATA_C1)

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phiC1(int p, double x, double* dphi);

/*
Compute the coefficients for the interpolant which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of interpolant
g2fg is (p + 1)-vector
*/
void compute_g2fgC1(short p, double* g2fg);

/*** DRIVER FUNCTIONS BELOW ***/

/*
*/
void driverC1(void);

/*
Tests phi and phi' without regard to method or scaling.
*/
void phi_test_allC1(void);

/*
Simple interface to output nesting coefficients to user.
*/
void print_nesting_coefficientsC1(void);

/*
*/
void phi_nesting_testC1(void);

// End of file