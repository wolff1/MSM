//-------|---------|---------|---------|---------|---------|---------|---------|
/*
centered b-spline phi
*/

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(int p, double x, double* dphi);

// End of file