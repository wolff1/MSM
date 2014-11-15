//-------|---------|---------|---------|---------|---------|---------|---------|
/*
b_spline.h - Structure for initializing, using, and uninitializing B-splines
*/

#ifndef	B_SPLINE_H
#define	B_SPLINE_H

#include "all.h"
#include "memory.h"
#include "interpolant.h"	//	remove?
#include "polynomial.h"

typedef struct
{
	//	Members
	INTERPOLANT		cmn;
	double*			omega;
	double*			omegap;
} B_SPLINE;

//	EXTERNAL Methods
void b_spline_initialize(void* Interpolant);
void b_spline_compute_g2g(void* Interpolant);
void b_spline_compute_tg2g(void* Interpolant);
void b_spline_evaluate(void* Interpolant);
void b_spline_uninitialize(void* Interpolant);

//	INTERNAL Methods
void b_spline_compute_omega(B_SPLINE* Bs);
void b_spline_compute_omega_prime(B_SPLINE* Bs);
void b_spline_compute_g2p(B_SPLINE* Bs);
void b_spline_compute_g2fg(B_SPLINE* Bs);

void compute_blurring_operator(short degree, double* B);
void compute_operator_inverse(short O_degree, const double* O,
							  short I_degree, double* I);
void convert_delta2_to_shifts(short cd_degree, const double* cd, double* s);
void compute_omega(short p, short n, double* omega);
void compute_omega_prime(short p, short mu, double* omegap);
double phi(short p, double x, double* dphi);
void new_phi(short p, double** g2p, short n, double* x, double* phi, double* dphi);
void compute_g2fg(short p, double* g2fg);

#endif

// End of file