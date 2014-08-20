//-------|---------|---------|---------|---------|---------|---------|---------|
/*
utility functions
*/

#define MAX(x,y) ((x > y) ? (x) : (y))
#define MIN(x,y) ((x < y) ? (x) : (y))

/*
Allocate and zero memory for 2D array of type double
*/
double** dynarr_d(int rows, int cols);

/*
Allocate and zero memory for vector (1D array) whose elements
have size "size"
*/
void* dynvec(int rows, size_t size);

/*
Free allocated memory
*/
void dynfree(void* ptr);

/*
Display elements of 2D array
*/
void display_dynarr_d(double** A, int rows, int cols);

/*
Display elements of vector
*/
void display_vector_d(double* x, int rows);

/*
a is input array of length a_degree+1 representing polynomial coefficients
b is input array of length b_degree+1 representing polynomial coefficients
c is output array of length a_degree+b_degree+1 representing polynomial
	coefficients s.t. c = a*b
NOTE:	It is assumed that the polynomials are in the same base
		By convention, use a for the higher degree, b for lower degree
*/
void mpoly(short a_degree, double* a, short b_degree, double* b, double* c);

/*
multiply polynomials where result has restricted degree

Same as above except that c is not to exceed c_degree.
NOTE: It is assumed that c has length c_degree+1
*/
void mpolyr(short a_degree, double* a, short b_degree, double* b, short c_degree, double* c);

// End of file