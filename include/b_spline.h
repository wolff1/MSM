//-------|---------|---------|---------|---------|---------|---------|---------|
/*
b_spline.h - Structure for initializing, using, and uninitializing B-splines
*/

#ifndef	B_SPLINE_H
#define	B_SPLINE_H

#include "all.h"
#include "memory.h"
#include "interpolant.h"
#include "polynomial.h"

typedef struct
{
	//	Members
	INTERPOLANT		cmn;
	short			mu;
	double*			omega;
	double*			omegap;
	STENCIL*		GammaIntermediate;
	STENCIL*		GammaTop;
} B_SPLINE;

//	EXTERNAL Methods
void b_spline_initialize(void* Interpolant, MSM_PARAMETERS* MsmParams);
void b_spline_compute_g2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams);
void b_spline_compute_tg2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams);
void b_spline_evaluate(void* Interpolant, long Len, double* X, double* F, double* DF);
void b_spline_uninitialize(void* Interpolant);

//	INTERNAL Methods
void b_spline_compute_omega(B_SPLINE* Bs);
void b_spline_compute_omega_prime(B_SPLINE* Bs);
void b_spline_compute_g2p(B_SPLINE* Bs);
void b_spline_compute_g2fg(B_SPLINE* Bs);

void b_spline_compute_blurring_operator(short degree, double* B);
void b_spline_compute_operator_inverse(short O_degree, const double* O, short I_degree, double* I);
void b_spline_convert_delta2_to_shifts(short cd_degree, const double* cd, double* s);

#endif

// End of file