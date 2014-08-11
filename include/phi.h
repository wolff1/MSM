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
Tests phi and phi' without regard to method or scaling.
*/
void phi_test_all(void);

// End of file