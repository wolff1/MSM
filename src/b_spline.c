//-------|---------|---------|---------|---------|---------|---------|---------|
/*
b_spline.c - Routines for initializing, using, and uninitializing B-splines
*/

#include "b_spline.h"

//	EXTERNAL Methods
void b_spline_initialize(void* Interpolant, MSM_PARAMETERS* MsmParams)
{
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;
	assert(Bs != NULL);
	assert(MsmParams != NULL);
//	printf("\tB_SPLINE initialization!\n");

	//	Set up COMMON members
	Bs->cmn.Size = sizeof(B_SPLINE);

	//	Set up COMMON function pointers
	Bs->cmn.copy = &b_spline_copy;
	Bs->cmn.compute_g2g = &b_spline_compute_g2g;
	Bs->cmn.compute_tg2g = &b_spline_compute_tg2g;
	Bs->cmn.evaluate = &b_spline_evaluate;
	Bs->cmn.uninitialize = &b_spline_uninitialize;

	//	Initialize Members
	Bs->mu = MsmParams->mu;
	Bs->omega = NULL;
	Bs->omegap = NULL;
	Bs->GammaIntermediate = NULL;
	Bs->GammaTop = NULL;

	//	Set up the B_SPLINE interpolant
	b_spline_compute_g2p(Bs);
	b_spline_compute_g2fg(Bs);
	b_spline_compute_omega(Bs);
	b_spline_compute_omega_prime(Bs);
//	b_spline_compute_g2g(Bs);		//	Happens in preprocess
//	b_spline_compute_tg2g(Bs);		//	Happens in preprocess
}

void b_spline_copy(void* Dst, void* Src)
{
	assert(Dst != NULL);
	assert(Src != NULL);

	//	--> INTERPOLANT is copied in interpolant_copy()

	//	Copy Members
	((B_SPLINE*)Dst)->mu = ((B_SPLINE*)Src)->mu;

	((B_SPLINE*)Dst)->omega = (double*) dynvec(((B_SPLINE*)Dst)->cmn.p/2+((B_SPLINE*)Dst)->mu+1, sizeof(double));
	memcpy(((B_SPLINE*)Dst)->omega, ((B_SPLINE*)Src)->omega, (((B_SPLINE*)Dst)->cmn.p/2+((B_SPLINE*)Dst)->mu+1)*sizeof(double));

	((B_SPLINE*)Dst)->omegap = (double*) dynvec(((B_SPLINE*)Dst)->cmn.p/2+((B_SPLINE*)Dst)->mu+1, sizeof(double));
	memcpy(((B_SPLINE*)Dst)->omegap, ((B_SPLINE*)Src)->omegap, (((B_SPLINE*)Dst)->cmn.p/2+((B_SPLINE*)Dst)->mu+1)*sizeof(double));

	//	The following two stencils do not survive the preprocessing step and don't need to be copied
//	((B_SPLINE*)Dst)->GammaIntermediate = (STENCIL*) dynmem(sizeof(STENCIL));
//	stencil_copy(((B_SPLINE*)Dst)->GammaIntermediate, ((B_SPLINE*)Src)->GammaIntermediate);

//	((B_SPLINE*)Dst)->GammaTop = (STENCIL*) dynmem(sizeof(STENCIL));
//	stencil_copy(((B_SPLINE*)Dst)->GammaTop, ((B_SPLINE*)Src)->GammaTop);
}

void b_spline_compute_g2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams)
{
	long			Size = 0;
	double			Scale = 1.0;
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;

	assert(Bs != NULL);
	assert(Softener != NULL);
	assert(MsmParams != NULL);
//	printf("\tB_SPLINE compute_g2g\n");

	//	Pre-processing (Intermediate levels)
	Size = (long) ceil(2.0*MsmParams->alpha);
//printf("Size = %ld\t", Size);
//Size = (long) ceil(1.5*Size);
//printf("1.5*Size = %ld\n", Size);
	Scale = MsmParams->h/MsmParams->a;

	//	Compute Gamma sequence (defined by theta function) for intermediate level grid(s)
	Bs->GammaIntermediate = (STENCIL*) dynmem(sizeof(STENCIL));
//	stencil_initialize(Bs->GammaIntermediate, Size, STENCIL_SHAPE_SPHERE);
	stencil_initialize(Bs->GammaIntermediate, Size, STENCIL_SHAPE_CUBE);
	stencil_populate(Bs->GammaIntermediate, Softener, STENCIL_FUNCTION_TYPE_THETA, Scale);
//	stencil_display(Bs->GammaIntermediate, MsmParams->h/MsmParams->a);

	//	Compute K^l by shifting Gamma sequence
	Bs->cmn.g2g = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_initialize(Bs->cmn.g2g, Bs->GammaIntermediate->Size, Bs->GammaIntermediate->Shape);
	stencil_shift(Bs->GammaIntermediate, Bs->cmn.p/2 + Bs->mu, Bs->omegap, Bs->cmn.g2g);

	stencil_free(Bs->GammaIntermediate);	//	We only need g2g, so its OK to delete Gamma sequence (FIXME - Save to file in case we need to enlarge later?)
	dynfree(Bs->GammaIntermediate);
//	stencil_display(Bs->cmn.g2g, 1.0);

//	stencil_naive(p, a, h, p_2+mu, omegap, k, c, &g2g);
}

void b_spline_compute_tg2g(void* Interpolant, SOFTENER* Softener, MSM_PARAMETERS* MsmParams)
{
	long			Size = 0;
	double			Scale = 1.0;
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;

	assert(Bs != NULL);
	assert(Softener != NULL);
	assert(MsmParams != NULL);
//	printf("\tB_SPLINE compute_tg2g\n");

	//	Pre-processing (Top level)
	Size = (long) 2*ceil(2.0*MsmParams->D/(pow(2.0,MsmParams->L)*MsmParams->h)) + MsmParams->p;	//	FIXME-Decouple from grid?
//printf("Size = %ld\t", Size);
//Size = (long) ceil(1.5*Size);
//printf("1.5*Size = %ld\n", Size);
	Scale = MsmParams->h/MsmParams->a;

	//	Compute Gamma sequence (defined by gamma function) for top level grid
	Bs->GammaTop = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_initialize(Bs->GammaTop, Size + (Bs->cmn.p/2 + Bs->mu), STENCIL_SHAPE_CUBE);
	stencil_populate(Bs->GammaTop, Softener, STENCIL_FUNCTION_TYPE_GAMMA, Scale);
//	stencil_display(Bs->GammaTop, MsmParams->h/MsmParams->a);

	//	Compute K^L by shifting Gamma sequence
	Bs->cmn.tg2g = (STENCIL*) dynmem(sizeof(STENCIL));
	stencil_initialize(Bs->cmn.tg2g, Size, Bs->GammaTop->Shape);
	stencil_shift_infinite(Bs->GammaTop, Bs->cmn.p/2 + Bs->mu, Bs->omegap, Bs->cmn.tg2g);

	stencil_free(Bs->GammaTop);	//	We only need tg2g, so its OK to delete Gamma sequence (FIXME - Save to file in case we need to enlarge later?)
	dynfree(Bs->GammaTop);
//	stencil_display(Bs->cmn.tg2g, 1.0);

//	stencil_naive_top(p, a, h, p_2+mu, omegap, k, c, &tg2g);
}

void b_spline_evaluate(void* Interpolant, long Len, double* X, double* F, double* DF)
{
    long		i = 0;
    short		j = 0;
	short		p = 0;
    short		p_2 = 0;
    short		k = 0;
    double		xp = 0.0;
    double		s = 0.0;
	double**	g2p = NULL;
	B_SPLINE*	Bs = (B_SPLINE*) Interpolant;

	assert(Bs != NULL);
    assert(X != NULL);
    assert(F != NULL);
    assert(DF != NULL);
//	printf("\tB_SPLINE evaluate\n");

	g2p = Bs->cmn.g2p;
	p = Bs->cmn.p;
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

void b_spline_uninitialize(void* Interpolant)
{
	B_SPLINE*		Bs = (B_SPLINE*) Interpolant;
	assert(Bs != NULL);
//	printf("\tUn-initializing B_SPLINE!\n");

	//	Free the dynamically allocated memory
	if (Bs->cmn.g2g)
		stencil_free(Bs->cmn.g2g);
	if (Bs->cmn.tg2g)
		stencil_free(Bs->cmn.tg2g);
	dynfree(Bs->cmn.g2g);
	dynfree(Bs->cmn.tg2g);
	dynfree(Bs->cmn.g2p[0]);
	dynfree(Bs->cmn.g2p);
	dynfree(Bs->cmn.g2fg);
	dynfree(Bs->omega);
	dynfree(Bs->omegap);
}

//	INTERNAL Methods
void b_spline_evaluate_recurrence(B_SPLINE* Bs, double x, double* f, double* df)
{
	short		p = 0;
	double**	Q = NULL;
	double**	Qp = NULL;
	double		u = 0.0;
	short		k = 0;
	short		v = 0;

	assert(Bs != NULL);
	p = Bs->cmn.p;

	// Initialize return variables
	*f = 0.0;
	*df = 0.0;

	if (fabs(x) < (double)(p/2.0))
	{
		// Dynamically create memory for Q and Qp 2D arrays
		Q = dynarr_d(p,p);		// Q(u)
		Qp = dynarr_d(p,p);		// Q'(u)
		u = x + (double)p/2.0;	// b/c centered B-spline

		// Start with k = 1 (Q_1 is the indicator function on [0,1[)
		k = 1;
		for (v = 0; v <= p-k; v++)
		{
			Q[0][v] = ((u-(double)v) < 1.0 && (u-(double)v) >= 0.0);
		}

		// use recurrance to build up to k=p-1
		for (k = 2; k <= p; k++)
		{
			for (v = 0; v <= p-k; v++)
			{
				Qp[k-1][v] = Q[k-2][v] - Q[k-2][v+1];
				Q [k-1][v] = (k*Q[k-2][v+1] + (u-(double)v)*Qp[k-1][v])/(k-1);
			}
		}

		// Set return variables
		*f = Q[p-1][0];
		*df = Qp[p-1][0];

		// Free dynamically allocated memory for Q and Qp 2D arrays
		dynfree(Q[0]);		// Free the data buffer memory
		dynfree(Qp[0]);
		dynfree(Q);			// Free the pointer-pointer memory
		dynfree(Qp);
	}
}

void b_spline_compute_g2p(B_SPLINE* Bs)
{
    double**    g2p = NULL;
    short       p = 0;
    short       p_2 = 0;
    short       i = 0;
    short       j = 0;
    short       k = 0;
    double      x = 0.0;
    double      xmk = 0.0;
    double      dx = 0;
    double**    A = NULL;
    double*     b = NULL;
	lapack_int  rc = 0;
	lapack_int* piv = NULL;

	assert(Bs != NULL);
//	printf("\tB_SPLINE compute_g2p\n");

	p = Bs->cmn.p;
	p_2 = p/2;
    dx = 1.0/(p_2-2+1);

	g2p = (double**) dynarr_d(p_2,p);
    A = (double**) dynarr_d(p,p);
    b = (double*) dynvec(p, sizeof(double));

    for (k = 0; k < p_2; k++)
    {
        x = k;
        for (j = 0; j < p_2; j++)
        {
            //  down column -> different values of x
//            printf("k=%hd, j=%hd, x=%f\n", k, j, x);
            xmk = x - (k+1);

            //  across powers in row
            A[j][0] = 1.0;
            A[j+p_2][0] = 0.0;
            A[j+p_2][1] = 1.0;
            for (i = 1; i < p-1; i++)
            {
                A[j][i] = A[j][i-1]*xmk;
                A[j+p_2][i+1] = (i+1)*A[j][i];
            }
            A[j][p-1] = A[j][p-2]*xmk;

            //  right hand side of linear system
//            b[j] = phi(p, x, &b[j+p_2]);
			b_spline_evaluate_recurrence(Bs, x, &b[j], &b[j+p_2]);
            x += dx;
        }

        //  Solve Ax = b with x being g2p[k]
        // NOTE: b is overwritten with x, A is overwritten with LU
        piv = (lapack_int*) dynvec(p,sizeof(lapack_int));
        rc = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) p, (lapack_int) 1,
                            A[0], (lapack_int) p, piv, b, (lapack_int) 1);
        dynfree(piv);
        assert(rc == 0); // zero is SUCCESS

        for (i = 0; i < p; i++)
        {
            g2p[k][i] = b[i];
        }
    }

	//	Point the B-spline object to the g2p coefficients
	Bs->cmn.g2p = g2p;
//	display_dynarr_d(Bs->cmn.g2p, Bs->cmn.p/2, Bs->cmn.p);

	//	Free dynamically allocated memory
	dynfree(A[0]);
    dynfree(A);
    dynfree(b);
}

void b_spline_compute_g2fg(B_SPLINE* Bs)
{
	double*		g2fg = NULL;
	short		n = 0;
	short		p = 0;
	short		p_2 = 0;
	long		num = 0;			// numerator
	long		den = 0;			// denominator

	assert(Bs != NULL);
//	printf("\tB_SPLINE compute_g2fg\n");

	p = Bs->cmn.p;
	p_2 = (short)p/2;
	num = p;
	den = p_2;

	g2fg = (double*) dynvec(p_2+1,sizeof(double));

	// Account for 2^(1-p) in denominator once
	for (n = 1; n < p; n++)
	{
		den *= 2;
	}

	// Build initial state of "p choose p/2"
	for (n = 1; n < p_2; n++) // exclude p_2 b/c it is in num and den
	{
		num *= (p-n);	// [(p)(p-1)(p-2)...(p/2+1)]
		den *= n;		// [(p/2)(p/2-1)...(2)(1)]*[2^(p-1)]
	}
	g2fg[0] = (double) num / (double) den;

	// Compute [p choose p/2 + n] with 1 <= n <= p/2
	for (n = 1; n <= p_2; n++)
	{
		g2fg[n] = (g2fg[n-1]*(double)(p_2-n+1)) / (double)(p_2+n);
	}

	//	Point B-spline object to newly computed nesting coefficients
	Bs->cmn.g2fg = g2fg;
//	display_vector_d(Bs->cmn.g2fg, Bs->cmn.p/2+1);
}

void b_spline_compute_omega(B_SPLINE* Bs)
{
	short		i = 0;
	short		p = 0;
	short		p_2 = 0;
	short		n = 10;
	double*		X = NULL;
	double*		B = NULL;
	double*		DB = NULL;
	double*		Psi = NULL;

	assert(Bs != NULL);
//	printf("\tB_SPLINE compute_omega\n");

	p = Bs->cmn.p;
	p_2 = p/2;
	n = p_2 + Bs->mu + 1;

	// Dynamically allocate memory
	X = (double*) dynvec(p_2, sizeof(double));
	B = (double*) dynvec(p_2, sizeof(double));
	DB = (double*) dynvec(p_2, sizeof(double));
	Psi = (double*) dynvec(1,sizeof(double));
	Bs->omega = (double*) dynvec(n, sizeof(double));

	//	Compute B(E)
	for (i = 0; i < p_2; i++)
	{
		X[i] = (double)i;
	}
	b_spline_evaluate(&Bs->cmn, p_2, X, B, DB);

	//	Set up Psi
	Psi[0] = 1.0;

	//	Solve for omega
    bibst_lss(-1, pow(2.0,-53), p_2, B, 1, 1, Psi, n, Bs->omega);
//	display_vector_d(Bs->omega, n);

	// Free dynamically allocated memory
	dynfree(X);
	dynfree(B);
	dynfree(DB);
	dynfree(Psi);
}

void b_spline_compute_omega_prime(B_SPLINE* Bs)
{
	short		i = 0;
	short		p = 0;
	short		p_2 = 0;
	short		mu = 0;
	double*		B = NULL;
	double*		B2 = NULL;
    double*     B2E = NULL;
	double*		A2 = NULL; // A^2 ~= B^{-2}
	double*		Ctmp = NULL;
	double*		C = NULL;
	double*		CE = NULL;
	double*		AE = NULL;
    double*     c = NULL;
	double*		pd = NULL;
	double*		pE = NULL;
	double*		c_full = NULL;
	double*		p_full = NULL;
	double*		high_terms = NULL;
	double*		omegap = NULL;

	assert(Bs != NULL);
//	printf("\tB_SPLINE compute_omega_prime\n");

	p = Bs->cmn.p;
	p_2 = p/2;
	mu = Bs->mu;

	// Dynamically allocate memory
	B = (double*) dynvec(p_2, sizeof(double));
	B2 = (double*) dynvec(p-1, sizeof(double));
    B2E = (double*) dynvec(p-1, sizeof(double));
	A2 = (double*) dynvec(p_2, sizeof(double));
	Ctmp = (double*) dynvec(p+p_2-2, sizeof(double));
	C = (double*) dynvec(p-2, sizeof(double));
	CE = (double*) dynvec(p-2, sizeof(double));
    c = (double*) dynvec(mu+1, sizeof(double));
	pd = (double*) dynvec(p_2+1,sizeof(double));
	pE = (double*) dynvec(p_2+1,sizeof(double));
	c_full = (double*) dynvec(2*mu+1,sizeof(double));
	p_full = (double*) dynvec(p+1,sizeof(double));
	high_terms = (double*) dynvec(p+2*mu+1,sizeof(double));
	omegap = (double*) dynvec(p_2+mu+1, sizeof(double));

	//	Compute B_p
	b_spline_compute_blurring_operator(p_2-1, B);

	//	Compute B^2
	mpoly(p_2-1, B, p_2-1, B, B2);

	//	Compute B^{-2} = 1 / B^2 = 1 / 1 - D => D = 1 - B^2
	b_spline_compute_operator_inverse(p-2, B2, p_2-1, A2);

	//	compute (A^2)(B^2) in order to get -C
	mpoly(p-2, B2, p_2-1, A2, Ctmp);
	//	keep only (delta^2)^{p/2} and higher
	for (i = p_2; i <= p+p_2-3; i++)
	{
		C[i-p_2] = -Ctmp[i];
	}

	//	convert C(delta^2) to C(E)
	b_spline_convert_delta2_to_shifts(p-3, C, CE);

	//	convert A^2(delta^2) to A(E) = omega'
	b_spline_convert_delta2_to_shifts(p_2-1, A2, omegap);

	// Solve bi-infinite linear system involving B^2(E) and C(E)
    b_spline_convert_delta2_to_shifts(p-2, B2, B2E);
    bibst_lss(-1, pow(2.0,-53), p-1, B2E, p-2, p-2, CE, mu+1, c);

	// use solution, c, to find (\delta^2)^{p/2} \sum c_m E^m in E basis
	pd[p_2] = 1.0;
	b_spline_convert_delta2_to_shifts(p_2,pd,pE);

	// convolve (c,pE) after forming full, symmetric operator
	c_full[mu] = c[0];
	for (i = 1; i <= mu; i++)
	{
		c_full[mu+i] = c[i];
		c_full[mu-i] = c[i];
	}

	p_full[p_2] = pE[0];
	for (i = 1; i <= p_2; i++)
	{
		p_full[p_2+i] = pE[i];
		p_full[p_2-i] = pE[i];
	}

	mpoly(2*mu, c_full, p, p_full, high_terms);

	//	add previous step to omega'
	for (i = 0; i <= p_2+mu; i++)
	{
		omegap[i] += high_terms[p_2+mu+i];
	}

	//	Set B-spline object to point to newly computed omega' values
	Bs->omegap = omegap;
//printf("\n");display_vector_d(Bs->omegap, p_2+mu+1);

	// Free dynamically allocated memory
	dynfree(B);
	dynfree(B2);
    dynfree(B2E);
	dynfree(A2);
	dynfree(Ctmp);
	dynfree(C);
	dynfree(CE);
    dynfree(c);
	dynfree(pd);
	dynfree(pE);
	dynfree(c_full);
	dynfree(p_full);
	dynfree(high_terms);
}

void b_spline_compute_blurring_operator(short degree, double* B)
{
/*
void b_spline_compute_blurring_operator(short degree)
where
	degree is the max degree of the operator in s^2 and
	B is the blurring operator.

To calculate blurring operator, B_p:
c_n = 1/a_0 (b_n - [a_1*c_{n-1} + a_2*c_{n-2} + ... + a_n*c_0])

NOTE:
	* we do not need to store b at all
	* factorials are not calculated anew in every iteration of the loops
		* f keeps a running product that is equivalent to factorial()
	* inner loop only touches columns which will change for given n,k
		* inner loop could be parallelized, although probably not worth it
	* for more details, go to
		http://pmsm.cs.purdue.edu/tiki-index.php?page=omegaprime
		or
		http://en.wikipedia.org/wiki/Formal_power_series#Dividing_series
*/
	double*		a = NULL;
	double**	c = NULL;
	double		f = 1.0;	// factorial
	short		n = 0;
	short		k = 0;
	short		i = 0;

	assert(degree > 0);
	assert(B != NULL);

	//	Allocate memory dynamically for computation
	a = (double*) dynvec(degree+1,sizeof(double));
	c = (double**) dynarr_d(degree,degree+1);

	// All but last iteration
	a[0] = 1.0;
	c[0][0] = 1.0;			// b_0 / a_0
	for (n = 1; n < degree; n++)
	{
		a[n] = -1.0 / (n*f);
		f = (2.0*n+1.0)*(2.0*n)*f;	// Keep running factorial
		c[n][n] = 1.0 / f;			// b_n
		for (k = 1; k <= n; k++)
		{
			for (i = 0; i <= n-k; i++)
			{
				c[n][i+k-1] = c[n][i+k-1] - a[k]*c[n-k][i];
			}
		}
	}

	// Last iteration, use B instead of c
	n = degree;
	a[n] = -1.0 / (n*f);
	f = (2.0*n+1.0)*(2.0*n)*f;	// Keep running factorial
	B[n] = 1.0 / f;				// c_n = b_n
	for (k = 1; k <= n; k++)
	{
		for (i = 0; i <= n-k; i++)
		{
			B[i+k-1] = B[i+k-1] - a[k]*c[n-k][i];
		}
	}

	//	Free dynamically allocated memory
	dynfree(a);
	dynfree(c[0]);
	dynfree(c);
}

void b_spline_compute_operator_inverse(short O_degree, const double* O, short I_degree, double* I)
{
/*
Use O^{-1} = 1 / O = 1 / 1 - D -> D = 1 - O
Then use geometric series
O^{-1} = 1 + D + D^2 + D^3 + ... + D^{I_degree}

NOTE:	O ("oh", not zero) is the input operator
		I ("eye") is the inverse, I = O^{-1}
*/
	short		i = 0;
	short		imax = MIN(O_degree,I_degree);
	short		k = 0;
	double*		D = NULL;
	double*		ID = NULL;

	assert(O_degree > 0);
	assert(O != NULL);
	assert(I_degree > 0);
	assert(I != NULL);

	// Dynamically allocate memory
	D = (double*) dynvec(O_degree+1, sizeof(double));
	ID = (double*) dynvec(I_degree+1, sizeof(double));

	// Initialize D = 1 - O
	D[0] = 0.0;	// 1 - O[0] where O[0] is 1
	for (i = 1; i <= O_degree; i++)
	{
		D[i] = -O[i];
	}

	// Ensure that memory for I is clean
	for (i = 0; i <= I_degree; i++)
	{
		I[i] = 0.0;
	}

	//	Set I using Horner's Rule
	//	First, I = 1 + D
	I[0] = 1.0;
	for (i = 1; i <= imax; i++)
	{
		I[i] = D[i];
	}

	// Then, iterate for each necessary power of D
	for (k = I_degree-1; k > 0; k--)
	{
		//I = I*D + 1;
		mpolyr(I_degree, I, O_degree, D, I_degree, ID);
		for (i = 0; i <= I_degree; i++)
		{
			I[i] = ID[i];
		}
		I[0] += 1.0;
	}

	// Free dynamically allocated memory
	dynfree(D);
	dynfree(ID);
}

void b_spline_convert_delta2_to_shifts(short cd_degree, const double* cd, double* s)
{
/*
Convert an operator from delta^2 to shifts: E-2+E

cd -> "central difference" operator (input)
s  -> "shifts" operator (output)

NOTE:
	- Both arrays should be length cd_degree+1
	- This routine does not initialize s
		- Instead it operates as s = s + ...

The kth term of (E-2+E)^n is given by:
s_k = ((-1)^(mod(n+k,2)))*nchoosek(2*n,k+n)
*/
	short		n = 0;
	short		k = 0;
	double		f = 1.0;	//	(2n)!
	double		nf = 1.0;	//	(n)!
	double		g = 1.0;	// coefficient for a term due to conversion
	double		sign1 = 1.0;
	double		sign2 = 1.0;

	assert(cd_degree > -1);
	assert(cd != NULL);
	assert(s != NULL);

	//	Convert cd(\delta^2) to s(E)
	for (n = 0; n <= cd_degree; n++)
	{
		g = f/(nf*nf);
		sign2 = sign1;
		s[0] = s[0] + sign2*cd[n]*g;
		for (k = 1; k <= n; k++)
		{
			sign2 = -sign2;
			g = ((n-k+1.0)*g) / ((double)(n+k));
			s[k] = s[k] + sign2*cd[n]*g;
		}

		// Update the (2n)! part of the combination for the next iteration
		f = (2.0*n+2.0)*(2.0*n+1.0)*f;
		nf = (n+1.0)*nf;
		sign1 = -sign1;
	}
}


// End of file