//-------|---------|---------|---------|---------|---------|---------|---------|
/*
memory.h - utility functions
*/

#include "all.h"

#define MEM_ALIGN	64

/*
Allocate and zero memory for 2D array of type double
*/
double** dynarr_d(unsigned long rows, unsigned long cols);

/*
Allocate and zero memory for vector (1D array) whose elements
have size "size"
*/
void* dynvec(unsigned long rows, size_t size);

/*
Free allocated memory
*/
void dynfree(void* ptr);

/*
Display elements of 2D array
*/
void display_dynarr_d(double** A, unsigned long rows, unsigned long cols);

/*
Display elements of vector
*/
void display_vector_d(double* x, unsigned long rows);

/*
bibst -> bi-infinite, banded, symmetric, Toeplitz
lss -> linear system solver

A_len   -> bandwidth
A       -> band (diagonal to edge of band)
b_len   -> length of rhs vector b, gives # of equations to solve
b_nnz   -> the number of non-zero elements in b
b       -> rhs vector
x_len   -> length of the solution vector x
x       -> solution vector (Ax = b)
*/
void bibst_lss(long max_itr, double tol,
				short A_len, double* A,
				short b_len, short b_nnz, double* b,
				short x_len, double* x);

// End of file