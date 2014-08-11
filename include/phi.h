//-------|---------|---------|---------|---------|---------|---------|---------|
/*
centered b-spline phi
*/

#define PHI_DATA		"phi.dat"

#define PHI_DATA_LEN	strlen(PHI_DATA)

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(int p, double x, double* dphi);

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of B-spline
g2fg is (p/2 + 1)-vector
*/
void compute_g2fg(short p, double* g2fg);

/*** DRIVER FUNCTIONS BELOW ***/

/*
Tests phi and phi' without regard to method or scaling.
*/
void phi_test_all(void);

/*
Simple interface to output nesting coefficients to user.
*/
void print_nesting_coefficients(void);

/*
*/
void phi_nesting_test(void);

// End of file