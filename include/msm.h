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
void msm_initialize(void* msm);
void msm_preprocess(void* msm);
void msm_evaluate(void* msm);
void msm_uninitialize(void* msm);

//	INTERNAL Methods
void short_range(MSM* msm);
void anterpolate(MSM* msm);
void restrict(MSM* msm);
void direct(MSM* msm);
void direct_top(MSM* msm);
void prolongate(MSM* msm);
void interpolate(MSM* msm);
void exclude(MSM* msm);

#endif

//	End of file