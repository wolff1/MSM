//-------|---------|---------|---------|---------|---------|---------|---------|
/*
msm.h - 
*/

#ifndef	MSM_H
#define	MSM_H

#include "all.h"
#include "interpolant.h"
#include "softener.h"
#include "method.h"

#include "even_powers.h"
#include "b_spline.h"
#include "c1_spline.h"

typedef struct
{
	double		a;
	double		h;
	short		p;
} MSM_PARAMETERS;

typedef struct
{
	char		ComputeShortRange;
	char		ComputeLongRange;
	char		ComputeExclusions;
	char		IsN;
	char		IsNLogN;
} MSM_OPTIONS;

typedef struct
{
	double**	Data;
} MSM_GRIDS;

typedef struct
{
	//	Members
	METHOD				cmn;
	MSM_PARAMETERS		prm;
	MSM_OPTIONS			opt;
	MSM_GRIDS			grd;
	INTERPOLANT			itp;
	SOFTENER			sft;
} MSM;

//	EXTERNAL Methods
void msm_initialize(void* Msm);
void msm_preprocess(void* Msm);
void msm_evaluate(void* Msm);
void msm_uninitialize(void* Msm);

//	INTERNAL Methods
void short_range(MSM* Msm);
void anterpolate(MSM* Msm);
void restrict(MSM* Msm);
void direct(MSM* Msm);
void direct_top(MSM* Msm);
void prolongate(MSM* Msm);
void interpolate(MSM* Msm);
void exclude(MSM* Msm);

#endif

//	End of file