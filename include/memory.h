//-------|---------|---------|---------|---------|---------|---------|---------|
/*
memory.h - utility functions
*/

#ifndef	MEMORY_H
#define	MEMORY_H

#include "all.h"

#define MEM_ALIGN	64

/*
Dynamically allocate memory and zero it out
*/
void* dynmem(size_t size);

/*
Free allocated memory
*/
void dynfree(void* ptr);

/*
Create 2D array of type double
*/
double** dynarr_d(unsigned long rows, unsigned long cols);

/*
Create vector (1D array) whose elements have size "size"
*/
void* dynvec(unsigned long rows, size_t size);

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

#endif

// End of file