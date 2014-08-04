//-------|---------|---------|---------|---------|---------|---------|---------|
/*
gamma.h - the softening function used to split the kernel
*/

#define GAMMA_DATA	"gamma.dat"
#define THETA_DATA	"theta.dat"

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
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
top_level is flag for intermediate or top level calculation
*/
double theta(double *c, int k, double x, double* dtheta, int long_range);

// End of file