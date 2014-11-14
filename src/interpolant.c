//-------|---------|---------|---------|---------|---------|---------|---------|
/*
centered b-spline phi
*/

#include "interpolant.h"

/*
	p gives the order of accuracy and p-1 is the degree of the interpolant
	x is the independent variable
	*dphi is output parameter which will contain the derivative of phi at x
*/
double phi(short p, double x, double* dphi)
{
	double f = 0.0;
	double df = 0.0;
	double** Q = NULL;
	double** Qp = NULL;
	double u = 0.0;
	short k = 0;
	short v = 0;

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
		f = Q[p-1][0];
		df = Qp[p-1][0];

		// Free dynamically allocated memory for Q and Qp 2D arrays
		dynfree(Q[0]);		// Free the data buffer memory
		dynfree(Qp[0]);
		dynfree(Q);			// Free the pointer-pointer memory
		dynfree(Qp);
	}

	if (dphi != NULL)
		*dphi = df;
	return f;
}

//  Horner's rule version of phi
void new_phi(short p, double** g2p, short n, double* x, double* phi, double* dphi)
{
    short       i = 0;
    short       j = 0;
    short       p_2 = p/2;
    short       k = 0;
    double      xp = 0.0;
    double      s = 0.0;

    assert(phi != NULL);
    assert(dphi != NULL);

    for (i = 0; i < n; i++)
    {
        phi[i] = 0.0;
        dphi[i] = 0.0;
        s = (x[i] > 0.0 ? 1.0 : -1.0);
        xp = fabs(x[i]);

        if (xp < p_2)
        {
            k = floor(xp);
            xp = xp - k - 1;

            phi[i] = g2p[k][p-1];
            dphi[i] = g2p[k][p-1]*(p-1);
            for (j = p-2; j > 0; j--)
            {
                phi[i] = phi[i]*xp + g2p[k][j];
                dphi[i] = dphi[i]*xp + g2p[k][j]*j;
            }
            phi[i] = phi[i]*xp + g2p[k][0];
            dphi[i] = dphi[i]*s;
        }
    }
}

/*
Compute the coefficients for the B-spines which allow them to be
nested from a grid to a finer grid.

p such that p-1 is degree of B-spline
g2fg is (p/2 + 1)-vector
*/
void compute_g2fg(short p, double* g2fg)
{
	short		n = 0;
	short		p_2 = (short)p/2;
	long		num = p;			// numerator
	long		den = p_2;			// denominator

	assert(g2fg != NULL);

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
}

//	End of file