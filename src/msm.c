//-------|---------|---------|---------|---------|---------|---------|---------|
/*
msm.c - 
*/

#include "msm.h"

//	EXTERNAL Methods
void msm_initialize(void* Method)
{
	size_t		Size = 0;
	void*		Init = NULL;
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	printf("Initializing MSM!\n");

	//	Set up COMMON function pointers
	Msm->cmn.preprocess = &msm_preprocess;
	Msm->cmn.evaluate = &msm_evaluate;
	Msm->cmn.uninitialize = &msm_uninitialize;

	//	Set up MSM parameters

	//	Set up MSM options

	//	Compute components which will not change throughout the simulation
/*	//	Initialize Interpolant
	if (1)
	{	//	B_SPLINE interpolant
		Size = sizeof(B_SPLINE);
		Init = &b_spline_initialize;
	}
	interpolant_initialize(&Msm->itp, Size, Init);
*/
	//	Initialize Softener
	if (1)
	{	//	EVEN_POWERS softening (aka "Taylor")
		Size = sizeof(EVEN_POWERS);
		Init = &even_powers_initialize;
	}
	softener_initialize(&Msm->sft, Size, Init);
}

void msm_preprocess(void* Method)
{
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	printf("MSM Preprocessing!\n");

	//	Compute Interpolant coefficients
}

void msm_evaluate(void* Method)
{
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	printf("MSM Evaluation!\n");

	if (Msm->opt.ComputeShortRange)
	{
		short_range(Msm);
	}

	if (Msm->opt.ComputeLongRange)
	{
		if (Msm->opt.IsN)
		{
			anterpolate(Msm);
			restrict(Msm);
			direct(Msm);
			direct_top(Msm);
			prolongate(Msm);
			interpolate(Msm);
		}

		if (Msm->opt.IsNLogN)
		{
			anterpolate(Msm);
			direct_top(Msm);
			interpolate(Msm);
		}

		if (Msm->opt.ComputeExclusions)
		{
			exclude(Msm);
		}
	}
}

void msm_uninitialize(void* Method)
{
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	printf("Un-initializing MSM!\n");

	dynfree(Msm);
}

//	INTERNAL Methods
void short_range(MSM* msm)
{
}

void anterpolate(MSM* msm)
{
}

void restrict(MSM* msm)
{
}

void direct(MSM* msm)
{
}

void direct_top(MSM* msm)
{
}

void prolongate(MSM* msm)
{
}

void interpolate(MSM* msm)
{
}

void exclude(MSM* msm)
{
}

#if 0
void msm_init(void)
{
	//	input: accuracy parameters:	a, h, p, mu, k, stencil shape,
	//			other parameters:	L (domain size), lmax (number of grids)
	//			internal objects:	g2p, g2fg, p2p, omega', g2g, tg2g
	//			grids:				grid, 1 to lmax

	//	g2p (g2p is p/2 x p)
	//compute_g2p(short p, double* g2p);

	//	g2fg (len = p/2+1)
	//compute_g2fg(short p, double* g2fg);

	//	p2p (p2p is k+1)
	//gamma_init(short k, double* p2p);

	//	omega' (omegap is p/2+mu+1)
	//compute_omega_prime(short p, short mu, double* omegap);

	//	g2g
	//	tg2g
	msm_preprocess();
}

void msm_preprocess(void)
{
	//	input: omegap', alpha = a/h, L = domain size, K-top, p, p_2, mu

	//	if g2g == NULL
	//		compute intermediate stencil (only happens once)
	//stencil_initialize(&theta, (long) ceil(2.0*alpha), STENCIL_SHAPE_SPHERE);
	//stencil_populate(&theta, c, k, STENCIL_FUNCTION_TYPE_THETA, h/a);
	//stencil_initialize(&g2g, theta.size, theta.shape);
	//stencil_shift(&theta, p_2 + mu, omegap, &g2g);
	//	write theta to a file for future use and to free up some memory?

	//	if tg2g == NULL or domain enlarged
	//		first full computation is resizing from 0 to X
	//		additional computations are resizing from X to Y
	//stencil_initialize(&gamma, 10 * (long) ceil(2.0*alpha), STENCIL_SHAPE_CUBE);
	//stencil_populate(&gamma, c, k, STENCIL_FUNCTION_TYPE_GAMMA, h/a);
	//stencil_initialize(&tg2g, gamma.size, gamma.shape);
	//stencil_shift(&gamma, p_2 + mu, omegap, &tg2g);
	//	write gamma to a file for future use and to free up some memory?
}

void msm_eval(void)
{
	//	input: r, q, options, parameters

	//	short_range
	//	anterpolate
	//		direct
	//		restrict
	//	direct_top
	//		prolongate
	//	interpolate
	//	exclusions
}
#endif

//	End of file