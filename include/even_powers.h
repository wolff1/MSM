//-------|---------|---------|---------|---------|---------|---------|---------|
/*
even_powers.h - the softening function used to split the kernel
*/

#define GAMMA_DATA			"gamma.dat"
#define THETA_DATA			"theta.dat"

#define GAMMA_DATA_LEN		strlen(GAMMA_DATA)
#define THETA_DATA_LEN		strlen(THETA_DATA)

#define	STENCIL_STORAGE(L)		(L+1)*(L+2)*(L+3)/6
#define	STENCIL_MAP_Z(z)		(z)*(z+1)*(z+2)/6
#define	STENCIL_MAP_Y(y)		(y)*(y+1)/2
#define STENCIL_MAP_X(x)		(x)

//	The following are for the "stacked" symmetric plane stencils
//	which are intermediate results while applying the anti-blurring
//	operator.
#define	STENCIL_MAP_Z2(L,z)		(z)*(L+2)*(L+1)/2	//	L is K->size
#define	STENCIL_MAP_Y2(y)		(y)*(y+1)/2
#define STENCIL_MAP_X2(x)		(x)

#define	STENCIL_STORAGE_2D(L)	(L+1)*(L+2)/2	//	Used for loop ranges

#define	STENCIL_SHAPE_SPHERE	1
#define	STENCIL_SHAPE_CUBE		2

#define	STENCIL_FUNCTION_TYPE_THETA		1
#define	STENCIL_FUNCTION_TYPE_GAMMA		2

typedef struct
{
	short			shape;
	long			size;	//	= zmax
	long*			ymax;
	long*			xmax;
	double*			data;
} STENCIL;

/*
Solves for coefficients of gamma
k is degree of continuity
x is k+1 vector to hold coefficients
*/
void gamma_init(short k, double* x);

/*
Evaluate gamma and gamma' at position x
c is coefficient vector for gamma
*/
double gamma(double *c, short k, double x, double* dgamma);

/*
Evaluate theta* and theta*' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta_star(double *c, short k, double x, double* dtheta_star);

/*
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta(double *c, short k, double x, double* dtheta);

/*
Evaluate theta and theta' at position x as single polynomial
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double thetap(double *c, short k, double x, double* dtheta);

/*
stencil operations below:
*/
void stencil_initialize(STENCIL* s, long size, short shape);
void stencil_populate(STENCIL* s, double* c, short k, short function_type,
						double h_a);
void stencil_display(STENCIL* s, double h_a);
void stencil_shift(STENCIL* s, short degree, double* omega_prime, STENCIL* K);
void stencil_free(STENCIL* s);
void build_Gamma_stencil(short k, double* c, double a, double h, double alpha, STENCIL* Gamma);

/*** DRIVER FUNCTIONS BELOW ***/

/*
Compare theta and thetap
*/
void test_thetas(void);

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
a:			cut-off distance
d:			length of domain starting at 0.0
X:			independent variable, 0.0 <= X <= d (vector, 1 x samples)
F:			dependent variable(s) (array, nlev x samples)
*/
void plot_splittings(int samples, int nlev, double a, double d,
						double* X, double** F);

/*
Tests gamma, theta, and theta_star without regard to method or scaling
*/
void gamma_test_all(void);

/*
Builds Gamma, anti-blurring operator, etc
*/
void test_preprocessing(void);

// End of file