//-------|---------|---------|---------|---------|---------|---------|---------|
/*
gamma.h - the softening function used to split the kernel
*/

/*
Solves for coefficients of gamma
k is degree of continuity
c is k+1 vector to hold coefficients
*/
void gamma_init(int k, double* x);

/*
Evaluate gamma and gamma' at position x
c is coefficient vector for gamma
*/
double gamma(double *c, int k, double x, double* dgamma);

// End of file