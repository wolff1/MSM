//-------|---------|---------|---------|---------|---------|---------|---------|
/*
c1_spline.c - 
*/

#include "c1_spline.h"

//	EXTERNAL Methods
void c1_spline_initialize(void* Interpolant, MSM_PARAMETERS* MsmParams)
{
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;
	assert(C1 != NULL);
	assert(MsmParams != NULL);
//	printf("\tC1_SPLINE initialization!\n");

	//	Set up COMMON members
	C1->cmn.Size = sizeof(C1_SPLINE);

	//	Set up COMMON function pointers
	C1->cmn.copy = &c1_spline_copy;
	C1->cmn.compute_g2g = &c1_spline_compute_g2g;
	C1->cmn.compute_tg2g = &c1_spline_compute_tg2g;
	C1->cmn.evaluate = &c1_spline_evaluate;
	C1->cmn.uninitialize = &c1_spline_uninitialize;

	//	Set up the C1_SPLINE interpolant
	c1_spline_compute_g2p(C1);
	c1_spline_compute_g2fg(C1);
}

void c1_spline_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	--> INTERPOLANT is copied in interpolant_copy()
}

void c1_spline_compute_g2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams)
{
	long			Size = 0;
	double			Scale = 1.0;
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;

	assert(C1 != NULL);
	assert(Softener != NULL);
	assert(MsmParams != NULL);
//	printf("\tC1_SPLINE compute_g2g\n");

	//	Pre-processing (Intermediate levels)
	Size = (long) ceil(2.0*MsmParams->alpha);
	Scale = MsmParams->h/MsmParams->a;

	//	Compute K^l sequence (defined by theta function) for intermediate level grid(s)
	C1->cmn.g2g = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_initialize(C1->cmn.g2g, Size, STENCIL_SHAPE_SPHERE);
	stencil_populate(C1->cmn.g2g, Softener, STENCIL_FUNCTION_TYPE_THETA, Scale);
//	stencil_display(C1->cmn.g2g, MsmParams->h/MsmParams->a);
}

void c1_spline_compute_tg2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams)
{
	long			Size = 0;
	double			Scale = 1.0;
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;

	assert(C1 != NULL);
	assert(Softener != NULL);
	assert(MsmParams != NULL);
//	printf("\tC1_SPLINE compute_tg2g\n");

	//	Pre-processing (Top level)
	Size = (long) ceil(MsmParams->D);
	Scale = MsmParams->h/MsmParams->a;

	//	Compute K^L sequence (defined by gamma function) for top level grid
	C1->cmn.tg2g = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_initialize(C1->cmn.tg2g, Size, STENCIL_SHAPE_CUBE);
	stencil_populate(C1->cmn.tg2g, Softener, STENCIL_FUNCTION_TYPE_GAMMA, Scale);
//	stencil_display(C1->cmn.tg2g, MsmParams->h/MsmParams->a);
}

void c1_spline_evaluate(void* Interpolant, long Len, double* X, double* F, double* DF)
{
    long		i = 0;
    short		j = 0;
	short		p = 0;
    short		p_2 = 0;
    short		k = 0;
    double		xp = 0.0;
    double		s = 0.0;
	double**	g2p = NULL;
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;

	assert(C1 != NULL);
    assert(X != NULL);
    assert(F != NULL);
    assert(DF != NULL);
//	printf("\tC1_SPLINE evaluate\n");

	g2p = C1->cmn.g2p;
	p = C1->cmn.p;
	p_2 = p/2;

    for (i = 0; i < Len; i++)
    {
        F[i] = 0.0;
        DF[i] = 0.0;
        s = (X[i] > 0.0 ? 1.0 : -1.0);
        xp = fabs(X[i]);

        if (xp < p_2)
        {
            k = floor(xp);
            xp = xp - k - 1;

            F[i] = g2p[k][p-1];
            DF[i] = g2p[k][p-1]*(p-1);
            for (j = p-2; j > 0; j--)
            {
                F[i] = F[i]*xp + g2p[k][j];
                DF[i] = DF[i]*xp + g2p[k][j]*j;
            }
            F[i] = F[i]*xp + g2p[k][0];
            DF[i] = DF[i]*s;
        }
    }
}

void c1_spline_uninitialize(void* Interpolant)
{
	C1_SPLINE*		C1 = (C1_SPLINE*) Interpolant;
	assert(C1 != NULL);
//	printf("\tUn-initializing C1_SPLINE!\n");

	//	Free the dynamically allocated memory
	stencil_free(C1->cmn.g2g);
	stencil_free(C1->cmn.tg2g);
	dynfree(C1->cmn.g2g);
	dynfree(C1->cmn.tg2g);
	dynfree(C1->cmn.g2p[0]);
	dynfree(C1->cmn.g2p);
	dynfree(C1->cmn.g2fg);
}

//	INTERNAL Methods
void c1_spline_compute_g2p(C1_SPLINE* C1)
{
	double**		g2p = NULL;
    short			p = 0;
    short			p_2 = 0;

	assert(C1 != NULL);
//	printf("\tC1_SPLINE compute_g2p\n");

	p = C1->cmn.p;
	p_2 = p/2;
	g2p = (double**) dynarr_d(p_2,p);

	switch (C1->cmn.p)
	{
	case 4:		//	Cubic
		g2p[0][0] = 0.0 / 2.0;
		g2p[0][1] = -1.0 / 2.0;
		g2p[0][2] = 4.0 / 2.0;
		g2p[0][3] = 3.0 / 2.0;

		g2p[1][0] = 0.0 / -2.0;
		g2p[1][1] = 0.0 / -2.0;
		g2p[1][2] = 1.0 / -2.0;
		g2p[1][3] = 1.0 / -2.0;
		break;
	case 6:		//	Quintic
		g2p[0][0] = 0.0 / 12.0;
		g2p[0][1] = -8.0 / 12.0;
		g2p[0][2] = 18.0 / 12.0;
		g2p[0][3] = 7.0 / 12.0;
		g2p[0][4] = -12.0 / 12.0;
		g2p[0][5] = -5.0 / 12.0;

		g2p[1][0] = 0.0 / -24.0;
		g2p[1][1] = -2.0 / -24.0;
		g2p[1][2] = 11.0 / -24.0;
		g2p[1][3] = 7.0 / -24.0;
		g2p[1][4] = -11.0 / -24.0;
		g2p[1][5] = -5.0 / -24.0;

		g2p[2][0] = 0.0 / 24.0;
		g2p[2][1] = 0.0 / 24.0;
		g2p[2][2] = 2.0 / 24.0;
		g2p[2][3] = 1.0 / 24.0;
		g2p[2][4] = -2.0 / 24.0;
		g2p[2][5] = -1.0 / 24.0;
		break;
	case 8:		//	Septic
		g2p[0][0] = 0.0 / 144.0;
		g2p[0][1] = -108.0 / 144.0;
		g2p[0][2] = 192.0 / 144.0;
		g2p[0][3] = 67.0 / 144.0;
		g2p[0][4] = -144.0 / 144.0;
		g2p[0][5] = -38.0 / 144.0;
		g2p[0][6] = 24.0 / 144.0;
		g2p[0][7] = 7.0 / 144.0;

		g2p[1][0] = 0.0 / -240.0;
		g2p[1][1] = -36.0 / -240.0;
		g2p[1][2] = 102.0 / -240.0;
		g2p[1][3] = 68.0 / -240.0;
		g2p[1][4] = -125.0 / -240.0;
		g2p[1][5] = -39.0 / -240.0;
		g2p[1][6] = 23.0 / -240.0;
		g2p[1][7] = 7.0 / -240.0;

		g2p[2][0] = 0.0 / 720.0;
		g2p[2][1] = -12.0 / 720.0;
		g2p[2][2] = 88.0 / 720.0;
		g2p[2][3] = 43.0 / 720.0;
		g2p[2][4] = -110.0 / 720.0;
		g2p[2][5] = -38.0 / 720.0;
		g2p[2][6] = 22.0 / 720.0;
		g2p[2][7] = 7.0 / 720.0;

		g2p[3][0] = 0.0 / -720.0;
		g2p[3][1] = 0.0 / -720.0;
		g2p[3][2] = 12.0 / -720.0;
		g2p[3][3] = 4.0 / -720.0;
		g2p[3][4] = -15.0 / -720.0;
		g2p[3][5] = -5.0 / -720.0;
		g2p[3][6] = 3.0 / -720.0;
		g2p[3][7] = 1.0 / -720.0;
		break;
	case 10:	//	Nonic
		g2p[0][0] = 0.0 / 2880.0;
		g2p[0][1] = -2304.0 / 2880.0;
		g2p[0][2] = 3600.0 / 2880.0;
		g2p[0][3] = 1300.0 / 2880.0;
		g2p[0][4] = -2740.0 / 2880.0;
		g2p[0][5] = -557.0 / 2880.0;
		g2p[0][6] = 620.0 / 2880.0;
		g2p[0][7] = 130.0 / 2880.0;
		g2p[0][8] = -40.0 / 2880.0;
		g2p[0][9] =  -9.0 / 2880.0;

		g2p[1][0] = 0.0 / -1440.0;
		g2p[1][1] = -288.0 / -1440.0;
		g2p[1][2] = 576.0 / -1440.0;
		g2p[1][3] = 446.0 / -1440.0;
		g2p[1][4] = -757.0 / -1440.0;
		g2p[1][5] = -199.0 / -1440.0;
		g2p[1][6] = 194.0 / -1440.0;
		g2p[1][7] = 44.0 / -1440.0;
		g2p[1][8] = -13.0 / -1440.0;
		g2p[1][9] =  -3.0 / -1440.0;

		g2p[2][0] = 0.0 / 10080.0;
		g2p[2][1] = -384.0 / 10080.0;
		g2p[2][2] = 1424.0 / 10080.0;
		g2p[2][3] = 828.0 / 10080.0;
		g2p[2][4] = -1932.0 / 10080.0;
		g2p[2][5] = -567.0 / 10080.0;
		g2p[2][6] = 546.0 / 10080.0;
		g2p[2][7] = 132.0 / 10080.0;
		g2p[2][8] = -38.0 / 10080.0;
		g2p[2][9] =  -9.0 / 10080.0;

		g2p[3][0] = 0.0 / -40320.0;
		g2p[3][1] = -144.0 / -40320.0;
		g2p[3][2] = 1332.0 / -40320.0;
		g2p[3][3] = 520.0 / -40320.0;
		g2p[3][4] = -1813.0 / -40320.0;
		g2p[3][5] = -497.0 / -40320.0;
		g2p[3][6] = 518.0 / -40320.0;
		g2p[3][7] = 130.0 / -40320.0;
		g2p[3][8] = -37.0 / -40320.0;
		g2p[3][9] =  -9.0 / -40320.0;

		g2p[4][0] = 0.0 / 40320.0;
		g2p[4][1] = 0.0 / 40320.0;
		g2p[4][2] = 144.0 / 40320.0;
		g2p[4][3] = 36.0 / 40320.0;
		g2p[4][4] = -196.0 / 40320.0;
		g2p[4][5] = -49.0 / 40320.0;
		g2p[4][6] = 56.0 / 40320.0;
		g2p[4][7] = 14.0 / 40320.0;
		g2p[4][8] = -4.0 / 40320.0;
		g2p[4][9] =  -1.0 / 40320.0;
		break;
	default:
		printf("Order p = %d is not currently supported.\n", C1->cmn.p);
		assert(0 == 1);
	}

	//	Point the B-spline object to the g2p coefficients
	C1->cmn.g2p = g2p;
//	display_dynarr_d(C1->cmn.g2p, p_2, p);
}

void c1_spline_compute_g2fg(C1_SPLINE* C1)
{
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	short		i = 0;

	assert(C1 != NULL);
//	printf("\tC1_SPLINE compute_g2fg\n");

	X = (double*) dynvec(C1->cmn.p+1, sizeof(double));
	F = (double*) dynvec(C1->cmn.p+1, sizeof(double));
	DF = (double*) dynvec(C1->cmn.p+1, sizeof(double));

	for (i = 0; i <= C1->cmn.p; i++)
	{
		X[i] = (double)i/2.0;
	}
	c1_spline_evaluate((void*)&C1->cmn, C1->cmn.p+1, X, F, DF);

	//	Point C1 object to newly computed nesting coefficients
	C1->cmn.g2fg = F;
//	display_vector_d(C1->cmn.g2fg, C1->cmn.p+1);

	//	Free dynamically allocated memory
	dynfree(X);
	dynfree(DF);
}

//	End of file