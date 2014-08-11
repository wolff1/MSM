//-------|---------|---------|---------|---------|---------|---------|---------|
/*
gamma.h - the softening function used to split the kernel
*/

#define GAMMA_DATA			"gamma.dat"
#define THETA_DATA			"theta.dat"

#define GAMMA_DATA_LEN		strlen(GAMMA_DATA)
#define THETA_DATA_LEN		strlen(THETA_DATA)

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

/*
Evaluate theta* and theta*' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta_star(double *c, int k, double x, double* dtheta_star);

/*
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta(double *c, int k, double x, double* dtheta);

/*** DRIVER FUNCTIONS BELOW ***/

/*
Prompts user for samples, number of levels, degree of continuity for smoothing,
	cut-off value and domain length.
Then uses theta_star, theta, and gamma to compute (and plot) the smoothings for
	the method.
*/
void splitting_test(void);

/*
samples:	number of data points
nlev:		number of levels (3->2 grids)
k:			continuity of smoothing
a:			cut-off distance
d:			length of domain starting at 0.0
X:			independent variable, 0.0 <= X <= d (vector, 1 x samples)
F:			dependent variable(s) (array, nlev x samples)
*/
void plot_splittings(int samples, int nlev, int k, double a, double d,
						double* X, double** F);

/*
Tests gamma, theta, and theta_star without regard to method or scaling
*/
void gamma_test_all(void);

// End of file