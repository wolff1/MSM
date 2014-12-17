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
//#include "simulation_domain.h"
#include "particle.h"

#include "even_powers.h"
#include "b_spline.h"
#include "c1_spline.h"

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
	//	Members
	METHOD				cmn;
	MSM_PARAMETERS		prm;
	MSM_OPTIONS			opt;
	INTERPOLANT*		itp;
	SOFTENER*			sft;
} MSM;

//	EXTERNAL Methods
void msm_initialize(void* Msm);
void msm_copy(void* Dst, void* Src);
void msm_preprocess(void* Msm, double DomainRadius);
void msm_evaluate(void* Msm, SIMULATION_DOMAIN* Domain);
void msm_uninitialize(void* Msm);

//	INTERNAL Methods
void msm_short_range(MSM* Msm, SIMULATION_DOMAIN* Domain);
void msm_anterpolate(MSM* Msm, SIMULATION_DOMAIN* Domain, short Level);
void msm_restrict(MSM* Msm);
void msm_direct(MSM* Msm);
void msm_direct_top(MSM* Msm);
void msm_prolongate(MSM* Msm);
void msm_interpolate(MSM* Msm);
void msm_exclude(MSM* Msm);

//	INTERNAL HELPER Methods
void msm_short_range_bin_to_bin(MSM* Msm, SIMULATION_DOMAIN* Domain, long* Next, long Particle1, long Particle2, long MaxIterationCount);
void msm_short_range_compute_self(MSM* Msm, SIMULATION_DOMAIN* Domain, long* First, long* Next, long* ParticlesPerBin, long XBinCount, long YBinCount, long x, long y, long z);
void msm_short_range_compute_neighbor(MSM* Msm, SIMULATION_DOMAIN* Domain, long* First, long* Next, long* ParticlesPerBin, long XBinCount, long YBinCount, long x, long y, long z, long l, long m, long n);

//	INTERNAL TESTING Method Prototypes
void msm_short_range_naive(MSM* Msm, SIMULATION_DOMAIN* Domain);

#endif

//	End of file