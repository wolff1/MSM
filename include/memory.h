//-------|---------|---------|---------|---------|---------|---------|---------|
/*
memory.h - utility functions
*/

#ifndef	MEMORY_H
#define	MEMORY_H

#include "all.h"

#define MEM_ALIGN	64

void* dynmem(size_t size);
void dynfree(void* ptr);
void mkl_memory_check(void);
double** dynarr_d(unsigned long rows, unsigned long cols);
void* dynvec(unsigned long rows, size_t size);
void display_dynarr_d(double** A, unsigned long rows, unsigned long cols);
void display_vector_d(double* x, unsigned long rows);

//	FIXME - THIS DOES NOT BELONG HERE!
void bibst_lss(long max_itr, double tol,
				short A_len, double* A,
				short b_len, short b_nnz, double* b,
				short x_len, double* x);

#endif

// End of file