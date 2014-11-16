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
	void*		Ptr = NULL;
	MSM*		Msm = (MSM*) Method;

	assert(Msm != NULL);
	printf("Initializing MSM!\n");

	//	Initialize COMMON function pointers
	Msm->cmn.preprocess = &msm_preprocess;
	Msm->cmn.evaluate = &msm_evaluate;
	Msm->cmn.uninitialize = &msm_uninitialize;

	//	Initialize MSM parameters
	Msm->prm.a = 12.5;
	Msm->prm.h = 2.5;
	Msm->prm.alpha = Msm->prm.a / Msm->prm.h;
	Msm->prm.p = 4;
	Msm->prm.k = 4;
	Msm->prm.mu = 10;
	Msm->prm.D = 0.0;	//	Not known until preprocess/evaluate

	//	Initialize MSM options
	Msm->opt.ComputeExclusions = 1;
	Msm->opt.ComputeLongRange = 1;
	Msm->opt.ComputeShortRange = 1;
	Msm->opt.IsN = 1;
	Msm->opt.IsNLogN = 1;

	//	Initialize INTERPOLANT
	Ptr = NULL;
	Msm->itp = NULL;
	if (0)
	{	//	B_SPLINE interpolant
		Size = sizeof(B_SPLINE);
		Init = &b_spline_initialize;
	}
	else
	{	//	C1_SPLINE interpolant
		Size = sizeof(C1_SPLINE);
		Init = &c1_spline_initialize;
	}
	interpolant_initialize(&Ptr, Size, Init, &Msm->prm);
	Msm->itp = (INTERPOLANT*) Ptr;

	//	Initialize SOFTENER
	Ptr = NULL;
	Msm->sft = NULL;
	if (1)
	{	//	EVEN_POWERS softening (aka "Taylor")
		Size = sizeof(EVEN_POWERS);
		Init = &even_powers_initialize;
	}
	softener_initialize(&Ptr, Size, Init, Msm->prm.k);
	Msm->sft = (SOFTENER*) Ptr;
}

void msm_preprocess(void* Method)
{
	MSM*		Msm = (MSM*) Method;
	assert(Msm != NULL);
	printf("MSM Preprocessing!\n");

	//	Compute Interpolant coefficients
	Msm->prm.D = 10.0;
	(*Msm->itp->compute_g2g)(Msm->itp, Msm->sft, &Msm->prm);
	(*Msm->itp->compute_tg2g)(Msm->itp, Msm->sft, &Msm->prm);
}

void msm_evaluate(void* Method)
{
	MSM*		Msm = (MSM*) Method;
	assert(Msm != NULL);
	printf("MSM Evaluation!\n");

	if (Msm->opt.ComputeShortRange)
	{
		msm_short_range(Msm);
	}

	if (Msm->opt.ComputeLongRange)
	{
		if (Msm->opt.IsN)
		{
			msm_anterpolate(Msm);
			msm_restrict(Msm);
			msm_direct(Msm);
			msm_direct_top(Msm);
			msm_prolongate(Msm);
			msm_interpolate(Msm);
		}

		if (Msm->opt.IsNLogN)
		{
			msm_anterpolate(Msm);
			msm_direct_top(Msm);
			msm_interpolate(Msm);
		}

		if (Msm->opt.ComputeExclusions)
		{
			msm_exclude(Msm);
		}
	}
}

void msm_uninitialize(void* Method)
{
	MSM*		Msm = (MSM*) Method;
	assert(Msm != NULL);
	printf("Un-initializing MSM!\n");

	//	Uninitialize SOFTENER
	(*Msm->sft->uninitialize)(Msm->sft);

	//	Uninitialize INTERPOLANT
	(*Msm->itp->uninitialize)(Msm->itp);

	//	Free dynamic memory allocated for MSM object
	dynfree(Msm);
}

//	INTERNAL Methods
void msm_short_range(MSM* Msm)
{
}

void msm_anterpolate(MSM* Msm)
{
}

void msm_restrict(MSM* Msm)
{
}

void msm_direct(MSM* Msm)
{
}

void msm_direct_top(MSM* Msm)
{
}

void msm_prolongate(MSM* Msm)
{
}

void msm_interpolate(MSM* Msm)
{
}

void msm_exclude(MSM* Msm)
{
}

//	End of file