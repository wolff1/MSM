//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.h - C1 phi
*/

#ifndef	C1_SPLINE_H
#define	C1_SPLINE_H

#include "all.h"
#include "memory.h"
#include "interpolant.h"

typedef struct
{
	//	Members
	INTERPOLANT		cmn;
} C1_SPLINE;

//	EXTERNAL Methods
void c1_spline_initialize(void* Interpolant);
void c1_spline_compute_g2g(void* Interpolant);
void c1_spline_compute_tg2g(void* Interpolant);
void c1_spline_evaluate(void* Interpolant);
void c1_spline_uninitialize(void* Interpolant);

//	INTERNAL Methods
void c1_spline_compute_g2p(C1_SPLINE* C1);
void c1_spline_compute_g2fg(C1_SPLINE* C1);

double phiC1(int p, double x, double* dphi);
void compute_g2fgC1(short p, double* g2fg);

#endif

// End of file