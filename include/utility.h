//-------|---------|---------|---------|---------|---------|---------|---------|
/*
utility functions
*/

#define MAX(x,y) ((x > y) ? (x) : (y))

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

// End of file