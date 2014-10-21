//-------|---------|---------|---------|---------|---------|---------|---------|
/*
gamma.cpp - the softening function used to split the kernel
*/

#include "stdafx.h"

/*
Solves for coefficients of gamma
k is degree of continuity
b will return k+1 vector of coefficients
*/
void gamma_init(short k, double* x)
{
	double** A = NULL;	// Coefficient matrix
	double* b = x;		// Alias output parameter x
	int i = 0;
	int j = 0;
	lapack_int rc = 0;
	lapack_int* piv = NULL;

	// b will first act as RHS, then return coefficients to caller
	assert(b != NULL);

	// Allocate memory for A
	A = dynarr_d(k+1, k+1);

	// Build row zero of A and b
	b[0] = 1.0;
	for (i = 0; i <= k; i++)
	{
		A[0][i] = 1.0;
	}

	// Build rows 1 through k of A and b
	for (i = 1; i <= k; i++)
	{
		b[i] = -i*b[i-1];
		// cols (starting with diagonal or subdiagonal)
		for (j = ceil(i/2.0); j <= k; j++)
		{
			// NOTE: 2j-i+1 is the power of the term of one less derivative
			A[i][j] = (2*j-i+1)*A[i-1][j];
		}
	}
/*
	// Display A
	printf("A = \n");
	display_dynarr_d(A,k+1, k+1);

	// Display b
	printf("b = \n");
	display_vector_d(b, k+1);
*/
	// Solve system Ax = b
	// NOTE: b is overwritten with x, A is overwritten with LU
	piv = (lapack_int*) dynvec(k+1,sizeof(lapack_int));
	rc = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int) k+1, (lapack_int) 1,
						A[0], (lapack_int) k+1, piv, b, (lapack_int) 1);
	dynfree(piv);
	assert(rc == 0); // zero is SUCCESS
/*
	// Display c
	printf("c = \n");
	display_vector_d(b, k+1);
*/
	// Free allocated memory
	dynfree(A[0]);
	dynfree(A);
}

/*
Evaluate gamma and gamma' at (positive) position x
c is coefficient vector for gamma
*/
double gamma(double *c, short k, double x, double* dgamma)
{
	double f = 0.0;
	double df = 0.0;
	double xx = x*x;
	int i = 0;

	// gamma is a spline which becomes f(x) = 1/x for x >= 1.0
	if (x >= 1.0)
	{
		f = 1.0/x;
		df = -f*f;	//	-f*f = -1.0/(x*x);
	}
	else
	{
		// Use Horner's rule to evaluate polynomial
		f = c[k];			// Even powers
		df = c[k]*(2*k);	// Odd powers
		for (i = k-1; i >= 1; i--)
		{
			f = f*xx + c[i];
			df = df*xx + c[i]*(2*i);
		}
		f = f*xx + c[0];
		df = df*x;
	}

	if (dgamma != NULL)
		*dgamma = df;
	return f;
}

/*
Evaluate theta* and theta*' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta_star(double *c, short k, double x, double* dtheta_star)
{
    double f = 0.0;
    double df = 0.0;

    f = 1.0/x - gamma(c, k, x, &df);
    df = -1.0/(x*x) - df;

    if (dtheta_star != NULL)
        *dtheta_star = df;
    return f;
}

/*
Evaluate theta and theta' at position x
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double theta(double *c, short k, double x, double* dtheta)
{
	double f = 0.0;
	double df = 0.0;
	double df2 = 0.0;

    f = gamma(c, k, x, &df) - 0.5*gamma(c, k, 0.5*x, &df2);
    df = df - 0.25*df2;

	if (dtheta != NULL)
		*dtheta = df;
	return f;
}

/*
Evaluate theta and theta' at position x as single polynomial
c is coefficient vector for gamma
k is degree of continuity of gamma
*/
double thetap(double *c, short k, double x, double* dtheta)
{
	double	f = 0.0;
	double	df = 0.0;
	double	exp2 = 0.0;
	double	xx = x*x;
	short	i = 0;

	if (x <= 1.0)
	{
		exp2 = pow(2.0,2*k+1);
		// Use Horner's rule to evaluate polynomial
		f = ((exp2-1.0)/exp2)*c[k];			// Even powers
		df = ((exp2-1.0)/exp2)*c[k]*(2*k);	// Odd powers
		for (i = k-1; i >= 1; i--)
		{
			exp2 *= 0.25;//exp2 = pow(2.0,2*i+1);
			f = f*xx + ((exp2-1.0)/exp2)*c[i];
			df = df*xx + ((exp2-1.0)/exp2)*c[i]*(2*i);
		}
		f = f*xx + (0.5)*c[0];
		df = df*x;
	}
	else if (x < 2.0)
	{
		f = 1.0/x - 0.5*gamma(c,k,0.5*x,&df);
		df = -1.0/xx - 0.25*df;
	}

	if (dtheta != NULL)
		*dtheta = df;
	return f;
}

//-------|---------|---------|---------|---------|---------|---------|---------|
void stencil_initialize(STENCIL* s, long size, short shape)
{
	double			r = 0.0;
	double			rr = 0.0;
	long			z = 0;
	long			y = 0;
	long			zz = 0;
	long			yyzz = 0;
	long			zi_2d = 0;

	//	Initialize stencil
	s->shape = shape;
	s->size = size;
	s->data = (double*)
				dynvec(STENCIL_STORAGE(s->size), sizeof(double));
	s->ymax = (long*)	//	One max per z
				dynvec(s->size+1, sizeof(long));
	s->xmax = (long*)	//	One max per (y,z)
				dynvec(STENCIL_STORAGE_2D(s->size), sizeof(long));

	//	Default all ymax and xmax values to -1
	for (z = 0; z <= s->size; z++)
	{
		zi_2d = STENCIL_MAP_Y(z);
		s->ymax[z] = -1;
		for (y = 0; y <= s->ymax[z]; y++)
		{
			s->xmax[STENCIL_MAP_X(y) + zi_2d] = -1;
		}
	}

	//	Set up loop ranges
	if (s->shape == STENCIL_SHAPE_SPHERE)
	{
		//	SPHERIC
		rr = (double) s->size * s->size;
		for (z = 0; z < s->size; z++)	//	NOTE: Not going to s->size
		{
			zz = z*z;
			zi_2d = STENCIL_MAP_Y(z);
			s->ymax[z] = MIN(z, (long) floor(sqrt(rr - zz)));
			for (y = 0; y <= s->ymax[z]; y++)
			{
				yyzz = y*y + zz;
				s->xmax[STENCIL_MAP_X(y) + zi_2d] =
					MIN(y, (long) floor(sqrt(rr - yyzz)));
			}
		}
	}
	else
	{
		//	CUBIC
		for (z = 0; z <= s->size; z++)
		{
			zi_2d = STENCIL_MAP_Y(z);
			s->ymax[z] = z;
			for (y = 0; y <= s->ymax[z]; y++)
			{
				s->xmax[STENCIL_MAP_X(y) + zi_2d] = y;
			}
		}
	}
}

//-------|---------|---------|---------|---------|---------|---------|---------|
void stencil_populate(STENCIL* s, double* c, short k, short type, double h_a)
{
	double			(*f)(double*,short,double,double*);
	long			z = 0;
	long			y = 0;
	long			x = 0;
	long			zz = 0;
	long			zi = 0;
	long			zi_2d = 0;
	long			yyzz = 0;
	long			yi = 0;

	if (type == STENCIL_FUNCTION_TYPE_THETA)
		f = &theta;
	else
		f = &gamma;

	for (z = 0; z <= s->size; z++)
	{
		zz = z*z;
		zi = STENCIL_MAP_Z(z);
		zi_2d = STENCIL_MAP_Y(z);
		for (y = 0; y <= s->ymax[z]; y++)
		{
			yyzz = y*y + zz;
			yi = zi + STENCIL_MAP_Y(y);
			for (x = 0; x <= s->xmax[STENCIL_MAP_X(y) + zi_2d]; x++)
			{
				s->data[yi + STENCIL_MAP_X(x)] = (*f)(c, k, h_a*sqrt((double)x*x + yyzz), NULL);
			}
		}
	}
}

//-------|---------|---------|---------|---------|---------|---------|---------|
void stencil_display(STENCIL* s, double h_a)
{
	long				z = 0;
	long				y = 0;
	long				x = 0;
	long				idx = 0;
	double				d = 0.0;

	for (z = 0; z <= s->size; z++)
	{
		for (y = 0; y <= s->ymax[z]; y++)
		{
			for (x = 0; x <= s->xmax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
				d = h_a*sqrt((double) x*x + y*y + z*z);
				idx = STENCIL_MAP_X(x) + STENCIL_MAP_Y(y) + STENCIL_MAP_Z(z);
				printf("(%ld,%ld,%ld) -> %e -> %+e\n", x,y,z, d, s->data[idx]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void stencil_shift(STENCIL* s, short degree, double* omegap, STENCIL* K)
{
	long				z = 0;
	long				y = 0;
	long				x = 0;
	short				n = degree;
	long				xx = 0, yy = 0, zz = 0;
	long				m = 0;
	long				r = s->size;
	double*				tmp1 = NULL;
	double*				tmp2 = NULL;

	assert(s != NULL);
	assert(K != NULL);
	assert(s->size == K->size);

	//	Create memory for intermediate stencils
	tmp1 = (double*) dynvec((K->size+1)*(K->size+1)*(K->size+2)/2,sizeof(double));

	//	Apply anti-blurring operator to s in Z direction, i.e., (A_z)s
//    printf("Apply anti-blurring in Z direction!\n");
//***NOTE: COULD RESTRICT LOOPS TO SPHERIC INDEXES OF K, I THINK***
	for (z = 0; z <= K->size; z++)
	{
		zz = STENCIL_MAP_Z2(K->size,z);
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= y; x++)
			{
				xx = STENCIL_MAP_X2(x);
//				printf("(%02ld,%02ld,%02ld) -> %ld", x,y,z, zz+yy+xx);

				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(r-z,n); m++)
				{	//	(x,y,z+m) -> z+m <= r -> m <= r-z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(z-y,n); m++)
				{	//	(x,y,z-m) -> z-m >= y -> m <= z-y
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z-m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = z-y+1; m <= MIN(z-x,n); m++)
				{	//	(x,z-m,y) -> z-m >= x -> m <= z-x
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{	//	(z-m,x,y) -> z-m >= 0 -> m <= z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
//					printf(":%ld", m);
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{	//	(m-z,x,y) -> m-z <= x -> m <= x+z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
//					printf(":%ld", m);
				}
				for (m = z+x+1; m <= MIN(y+z,n); m++)
				{	//	(x,m-z,y) -> m-z <= y -> m <= y+z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{	//	(x,y,m-z) -> m-z <= r -> m <= z+r
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}
		}

		for (y = z+1; y <= K->size; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= z; x++)
			{
				xx = STENCIL_MAP_X2(x);
//***NOTE: x <= z < y
//				printf("(%02ld,%02ld,%02ld) -> %ld", x,y,z, zz+yy+xx);

				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z)+STENCIL_MAP_X(x)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = y-z+1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(z-x,n); m++)
				{	// z-m >= x -> m <= z-x
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{	//	z-m >= 0 ->  m <= z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
//					printf(":%ld", m);
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{	//	m-z <= x -> m <= x+z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
//					printf(":%ld", m);
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{	//	m-z <= y -> m <= y+z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{	//	m-z <= r -> m <= z+r
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}

			for (x = z+1; x <= y; x++)
			{
				xx = STENCIL_MAP_X2(x);
//***NOTE: z < x <= y
//				printf("(%02ld,%02ld,%02ld) -> %ld", x,y,z, zz+yy+xx);

				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(x-z,n); m++)
				{	// z+m <= x -> m <= x-z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z+m)];
//					printf(":%ld", m);
				}
				for (m = x-z+1; m <= MIN(y-z,n); m++)
				{	//	z+m <= y -> m <= y-z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = y-z+1; m <= MIN(r-z,n); m++)
				{	//	z+m <= r -> m <= r-z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(z,n); m++)
				{	// z-m >= 0 -> m <= z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
//					printf(":%ld", m);
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{	//	m-z <= x -> m <= x+z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
//					printf(":%ld", m);
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{	//	m-z <= y -> m <= y+z
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{	//	m-z <= r -> m <= z+r
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}
		}
	}
/*
	//	Display "stacked" half-plane where x-y plane is symmetric
	printf("tmp1:\n");
	for (z = 0; z <= K->size; z++)
	{
		zz = STENCIL_MAP_Z2(K->size,z);//z*(K->size+2)*(K->size+1)/2;
		for (y = 0; y <= K->size; y++)
		{
			yy = STENCIL_MAP_Y2(y);//y*(y+1)/2;
			for (x = 0; x <= y; x++)
			{
				xx = STENCIL_MAP_X2(x);//x;
				printf("%+e ", tmp1[zz+yy+xx]);
			}
			printf("\n");
		}
		printf("\n");
	}
*/
//    printf("Apply anti-blurring in Y direction!\n");
	tmp2 = (double*) dynvec((K->size+1)*(K->size+1)*(K->size+2)/2,sizeof(double));
	r = K->size;	//	This represents the length of the "stacked" symmetric plane stencils
	//	Apply anti-blurring operator to (A_z)s in Y direction, i.e., (A_y)(A_z)s
	for (z = 0; z <= K->size; z++)
	{
		zz = STENCIL_MAP_Y2(z);
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_X2(y);
			for (x = 0; x <= y; x++)
			{
//	x <= y <= z
				xx = STENCIL_MAP_Z2(K->size,x);
//				printf("(%02ld,%02ld,%02ld) -> %02ld", x,y,z, xx+zz+yy);

				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y)+STENCIL_MAP_X2(x)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(y-x,n); m++)
				{	//	y-m >= x -> m <= y-x
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y-m)+STENCIL_MAP_X2(x)];
//					printf(":%ld", m);
				}
				for (m = y-x+1; m <= MIN(y,n); m++)
				{	//	y-m >= 0 -> m <= y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
//					printf(":%ld", m);
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{	//	m-y <= x -> m <= x+y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
//					printf(":%ld", m);
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{	//	m-y <= r -> m <= r+y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}

			for (x = y+1; x <= z; x++)
			{
//	y <= x <= z
				xx = STENCIL_MAP_Z2(K->size,x);
//				printf("(%02ld,%02ld,%02ld) -> %02ld", x,y,z, xx+zz+yy);

				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(x-y,n); m++)
				{	//	y+m <= x -> m <= x-y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
//					printf(":%ld", m);
				}
				for (m = x-y+1; m <= MIN(r-y,n); m++)
				{	//	y+m <= r -> m <= r-y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(y,n); m++)
				{	//	y-m >= 0 -> m <= y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
//					printf(":%ld", m);
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{	//	m-y <= x -> m <= x+y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
//					printf(":%ld", m);
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{	//	m-y <= r -> m <= r+y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}

			for (x = z+1; x <= K->size; x++)
			{
//	y <= z <= x
				xx = STENCIL_MAP_Z2(K->size,x);
//				printf("(%02ld,%02ld,%02ld) -> %02ld", x,y,z, xx+zz+yy);

				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(z-y,n); m++)
				{	//	y+m <= z -> m <= z-y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y+m)];
//					printf(":%ld", m);
				}
				for (m = z-y+1; m <= MIN(r-y,n); m++)
				{	//	y+m <= r -> m <= r-y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(z)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(y,n); m++)
				{	//	y-m >= 0 -> m <= y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y-m)];
//					printf(":%ld", m);
				}
				for (m = y+1; m <= MIN(z+y,n); m++)
				{	//	m-y <= z -> m <= y+z
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(m-y)];
//					printf(":%ld", m);
				}
				for (m = y+z+1; m <= MIN(r+y,n); m++)
				{	//	m-y <= r -> m <= r+y
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(z)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}
		}
	}

	//	Free dynamically allocated memory for tmp1
	dynfree(tmp1);
/*
	//	Display "stacked" half-plane where y-z plane is symmetric
	printf("\ntmp2\n");
	for (x = 0; x <= K->size; x++)
	{
		//xx = x*(K->size+2)*(K->size+1)/2;
		xx = STENCIL_MAP_Z2(K->size,x);
		for (z = 0; z <= K->size; z++)
		{
			//zz = z*(z+1)/2;
			zz = STENCIL_MAP_Y2(z);
			for (y = 0; y <= z; y++)
			{
				//yy = y;
				yy = STENCIL_MAP_X2(y);
				printf("%+e ", tmp2[xx+zz+yy]);
			}
			printf("\n");
		}
		printf("\n");
	}
*/
//    printf("Apply anti-blurring in X direction!\n");
	r = K->size;
	//	Apply anti-blurring operator to (A_y)(A_z)s in X direction, i.e., (A_x)(A_y)(A_z)s
	for (z = 0; z <= K->size; z++)
	{
		zz = STENCIL_MAP_Z(z);
		for (y = 0; y <= K->ymax[z]; y++)
		{
			yy = STENCIL_MAP_Y(y);
			for (x = 0; x <= K->xmax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
				xx = STENCIL_MAP_X(x);
//				printf("(%02ld,%02ld,%02ld) -> %02ld", x,y,z, zz+yy+xx);

				K->data[zz+yy+xx] = omegap[0]*tmp2[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
//				printf(":%ld", 0);
				for (m = 1; m <= MIN(r-x,n); m++)
				{	//	x+m <= r -> m <= r-x
					K->data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->size,x+m)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
//					printf(":%ld", m);
				}
				for (m = 1; m <= MIN(x,n); m++)
				{	//	x-m >= 0 -> m <= x
					K->data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->size,x-m)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
//					printf(":%ld", m);
				}
				for (m = x+1; m <= MIN(r+x,n); m++)
				{	//	m-x <= r -> m <= r+x
					K->data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->size,m-x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
//					printf(":%ld", m);
				}
//				printf("\n");
			}
		}
	}

	//	Free dynamically allocated memory for tmp2
	dynfree(tmp2);
}

//-------|---------|---------|---------|---------|---------|---------|---------|
void stencil_free(STENCIL* s)
{
	if (s->data != NULL)
		dynfree(s->data);

	if (s->xmax != NULL)
		dynfree(s->xmax);

	if (s->ymax != NULL)
		dynfree(s->ymax);
}

/*** DRIVER FUNCTIONS BELOW ***/

/*
Compare theta and thetap
*/
void test_thetas(void)
{
	int		i = 0;
	int		samples = 1000;
	double	t1 = 0.0;
	double	dt1 = 0.0;
	double	t2 = 0.0;
	double	dt2 = 0.0;
	double	e = 0.0;
	double	de = 0.0;
	short	k = 0;
	double*	c = NULL;
	double	x = 0.0;

	printf("Enter samples: ");
	scanf("%d", &samples);

	printf("Enter k: ");
	scanf("%hd", &k);

	c = (double*) dynvec(k+1,sizeof(double));

	// Get gamma coefficients
	gamma_init(k, c);

	for (i = 0; i <= samples; i++)
	{
		x = ((double)i/(double)samples)*2.5;
		t1 = theta(c,k,x,&dt1);
		t2 = thetap(c,k,x,&dt2);

		e = fabs(t1 - t2);
		de = fabs(dt1 - dt2);

		printf("%f: Error = %e, dError = %e\n", x, e, de);
	}

	dynfree(c);
}

/*
Prompts user for samples, number of levels, degree of continuity for smoothing,
	cut-off value and domain length.
Then uses theta_star, theta, and gamma to compute (and plot) the smoothings for
	the method.
*/
void splitting_test(void)
{
	// Static memory variables
	int			i = 0;
	int			samples = 0;
	short		k = 0;
	double		a = 0.0;
	double		one_over_a = 0.0;
	double		al = 0.0;
	int			nlev = 0;
	int			l = 0;
	double		d = 0.0;
	double		f = 0.0;	//FIXME - Remove f, df, tol, and tol check below
	double		df = 0.0;
	double		tol = pow(0.1,15);
	// Dynamic memory variables
	double*		c = NULL;	// Gamma
	double*		X = NULL;
	double**	F = NULL;
	double**	DF = NULL;
	char*		data_file = NULL;

	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get nlev where nlev is the number of levels
	printf("Please enter levels parameter nlev = ");
	scanf("%d", &nlev);
	assert(nlev > 1);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%hd", &k);
	assert(k > 0);

	// Get a where a is cut-off distance
	printf("Please enter cut-off parameter a = ");
	scanf("%lf", &a);
	assert(a > 0.0);
	one_over_a = 1.0/a;

	// Get d where d is domain
	printf("Please enter domain parameter d = ");
	scanf("%lf", &d);
	assert(d > 0.0);

	// Dynamically allocate memory
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double**) dynarr_d(nlev+1,samples+1);
	DF = (double**) dynarr_d(nlev+1,samples+1);
	c = (double*) dynvec(k+1,sizeof(double));
	data_file = (char*) dynvec(GP_DATA_DIR_LEN +
						MAX(PHI_DATA_LEN,GAMMA_DATA_LEN) + 1, sizeof(char));
	assert(X != NULL);
	assert(F != NULL);
	assert(DF != NULL);
	assert(c != NULL);
	assert(data_file != NULL);

	// Initialize gamma coefficients
	gamma_init(k, c);

	// Evaluate different splitting functions over domain
	for (i = 0; i <= samples; i++)
	{
f = 0.0;	// sanity check
df = 0.0;	// sanity check
		X[i] = (d*(double)i/(double)samples);

		// Theta* - Short range part of splitting (finite)
		al = one_over_a;
		F[0][i] = al*theta_star(c, k, al*X[i], &DF[0][i]);
		DF[0][i] = al*al*DF[0][i];
f += F[0][i];
df += DF[0][i];

		for (l = 1; l < nlev - 1; l++)
		{
			// Theta - Intermediate long range part(s) of splitting (finite)
			F[l][i] = al*theta(c, k, al*X[i], &DF[l][i]);
			DF[l][i] = al*al*DF[l][i];
f += F[l][i];
df += DF[l][i];

			al = 0.5*al;
		}
		// Gamma - Top level long range part of splitting (infinit)
		F[l][i] = gamma(c, k, al*X[i], &DF[l][i]);
		F[l][i] = al*F[l][i];
		DF[l][i] = al*al*DF[l][i];
f += F[l][i];
df += DF[l][i];

		// Kernel: 1/X and Kernel' = -1/X^2
		F[nlev][i] = 1.0/X[i];
		DF[nlev][i] = -F[nlev][i]*F[nlev][i];
if (i > 0)
{
    double Frel = fabs(F[nlev][i] - f)/fabs(F[nlev][i]);
    double DFrel = fabs(DF[nlev][i] - df)/fabs(DF[nlev][i]);
	if (Frel >= tol)
		printf("i = %d, X = %f, F = %f, f = %f, |.| = %e\n",
				i, X[i], F[nlev][i], f, Frel);
	if (DFrel >= tol)
		printf("i = %d, X = %f, DF = %f, df = %f, |.| = %e\n",
				i, X[i], DF[nlev][i], df, DFrel);
	assert(Frel < tol);
	assert(DFrel < tol);
}
	}

	// Plot splittings on single graph along with f(x) = 1/x
	plot_splittings(samples, nlev, a, d, X, F);
//	plot_splittings(samples, nlev, k, a, d, X, DF);

	// Free dynamically allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F[0]);
	dynfree(F);
	dynfree(DF[0]);
	dynfree(DF);
	dynfree(data_file);
}

/*
samples:	number of data points
nlev:		number of levels (3->2 grids)
a:			cut-off distance
d:			length of domain starting at 0.0
X:			independent variable, 0.0 <= X <= d (vector, 1 x samples)
F:			dependent variable(s) (array, nlev x samples)
*/
void plot_splittings(int samples, int nlev, double a, double d,
						double* X, double** F)
{
	// Static memory variables
	char*	cmd_file_name = "splitting_";
	char*	cmd_file_extension = ".gp";
	size_t	cmd_file_name_len = strlen(cmd_file_name) + GP_TERM_LEN +
								strlen(cmd_file_extension);
	// File stuff
	size_t  bufmax = 0;
	size_t  buf2max = 0;
	size_t  buflen = 0;
	size_t  bytes = 0;
	int		i = 0;
	int		l = 0;
	// Dynamic memory variables
	FILE*	data = NULL;
	FILE*	cmd = NULL;
	char*	data_file = NULL;
	char*	cmd_file = NULL;
	char*	buf = NULL;
	char*	buf2 = NULL;
    char*   buf2_1 =    "set term %s\n"
                        "set xlabel 'r'\n"
						"set ylabel 'g_l(r)'\n"
						// below two lines rotates by -90 degrees
                        //"set lmargin 10\n"
                        //"set label 1 'g_l(r)' at graph -0.1, graph 0.5\n"
                        "set label 2 '(a)' at graph 0.39, graph -0.08\n"
                        "set label 3 '(2a)' at graph 0.785, graph -0.08\n"
                        "set title 'Kernel Splitting for %d-level MSM'\n"
                        "set grid\n"
//                      "set style data lines\n"
                        "set style data linespoints\n"
                        "set yrange [ 0.0 : %7.3f ]\n"
						"set key box\n"
                        "plot ";
    char*   buf2_2 =    "data_file using 1:%d title \"g_%d\" lc rgb \"black\",";
    char*   buf2_3 =    "data_file using 1:%d title \"1/r\" lc rgb \"black\"\n"
                        "pause -1\n"
                        "quit\n";

	assert(X != NULL);
	assert(F != NULL);

	// Dynamically allocate memory
	data_file = (char*) dynvec(GP_DATA_DIR_LEN + GP_DATA_TMP_LEN + 1,
								sizeof(char));
	cmd_file = (char*) dynvec(GP_CMD_DIR_LEN + cmd_file_name_len + 1,
								sizeof(char));

	// Build file name(s)
	sprintf(data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(cmd_file, "%s%s%s%s",
					GP_CMD_DIR, cmd_file_name, GP_TERM, cmd_file_extension);

	// Create DATA file with X, F[0], F[1], ..., F[nlev] columns and samples rows
	data = fopen(data_file, "w");
	assert(data != NULL);

	// Write DATA file
	bufmax = 64*(nlev+2);	// nlev + 2 columns, each a max of 64 chars wide 
	buf = (char*) dynvec(bufmax+1,sizeof(char));	// + 1 for NULL
	// Adjust F[0][0] and F[nlev][0] to be large number (Not necessary on OSX)
	F[0][0] = 1000.0;
	F[nlev][0] = 1000.0;
	for (i = 0; i <= samples; i++)
	{
		// Add X to buffer
		buflen = sprintf(buf, "%.32f", X[i]);
		assert(bufmax > buflen);

		for (l = 0; l < nlev; l++)
		{
			// Add F[l][i] to buffer
			buflen += sprintf(&buf[buflen], " %.32f", F[l][i]);
			assert(bufmax > buflen);
		}
		// Add F[nlev][i] to buffer
		buflen += sprintf(&buf[buflen], " %.32f\n", F[nlev][i]);
		assert(bufmax > buflen);

		// Write buffer to file
		bytes = fwrite(buf, sizeof(char), buflen, data);
		if (bytes < buflen)
		{
			printf("<%lu> bytes written to temp file <%s> instead of <%lu>\n",
				bytes, data_file, buflen);
			break;
		}

		// Clear out buffer for next usage
		memset(buf, 0, bufmax+1);
	}

	// Close DATA file (ensure buffer is flushed to disk)
	if (fclose(data))
	{
		printf("Error closing DATA file.\n");
	}

	// Create COMMAND file to plot all of the columns above
	cmd = fopen(cmd_file, "w");
	assert(cmd != NULL);

// Write CMD file
	buf2max = strlen(buf2_1) + nlev*strlen(buf2_2) + strlen(buf2_3) + 64;
	buf2 = (char*) dynvec(buf2max+1,sizeof(char));

	// Create COMMAND file buffer
	buflen = sprintf(buf2, buf2_1, GP_TERM, nlev, 1.33*F[1][0]);
	assert(buf2max > buflen);

	// Plot a line for each level (columns 1 -> nlev+1)
	for (i = 0; i < nlev; i++)
	{
		buflen += sprintf(&buf2[buflen], buf2_2, i+2, i);
		assert(buf2max > buflen);
	}

	// Plot line for 1/r (column nlev+2) and finish file
	buflen += sprintf(&buf2[buflen], buf2_3, nlev+2);
	assert(buf2max > buflen);

	// Write buffer to file
	bytes = fwrite(buf2, sizeof(char), buflen, cmd);
	if (bytes < buflen)
	{
		printf("<%lu> bytes written to temporary file <%s> instead of <%lu>\n",
				bytes, cmd_file, buflen);
	}

	// Close CMD file (ensure buffer is flushed to disk)
	if (fclose(cmd))
	{
		printf("Error closing DATA file.\n");
	}

	// Call gnuplot to plot DATA file using COMMAND file
	plotf2d(cmd_file, data_file);

	// Delete data file
	if (remove(data_file))
	{
		printf("Error removing DATA file <%s>.\n", data_file);
	}

	// Delete command data file
	if (remove(cmd_file))
	{
		printf("Error removing COMMAND file <%s>.\n", cmd_file);
	}

	// Free dynamically allocated memory
	dynfree(data_file);
	dynfree(cmd_file);
	dynfree(buf);
	dynfree(buf2);
}

/*
Tests gamma, theta, and theta_star without regard to method or scaling
*/
void gamma_test_all(void)
{
	int			i = 0;
	int			samples = 0;
	short		k = 0;
	double*		X = NULL;
	double*		F = NULL;
	double*		DF = NULL;
	double*		c = NULL;

	/**************************************************************************/
	// Get number of data points to record
	printf("Please enter the number of samples = ");
	scanf("%d", &samples);
	assert(samples > 0);

	// Get k where k is degree of continuity of the softener
	printf("Please enter softening parameter k = ");
	scanf("%hd", &k);
	assert(k > 0);

	// Create arrays for dependent and independent variables
	X = (double*) dynvec(samples+1, sizeof(double));
	F = (double*) dynvec(samples+1, sizeof(double));
	DF = (double*) dynvec(samples+1, sizeof(double));
	c = (double*) dynvec(k+1,sizeof(double));
	assert(X != NULL);
	assert(F != NULL);
	assert(DF != NULL);
	assert(c != NULL);

	/**************************************************************************/
	// Test softener
	gamma_init(k, c);
	for (i = 0; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = gamma(c, k, X[i], &DF[i]);
/*
		printf("%02d:\tgamma(%f) = %f\tgamma'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Gamma(x)
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, F);
	//	Show Gamma'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[0], X[samples]);
	plot2d(samples, X, DF);

	/**************************************************************************/
	// Test Theta* - short range splitting function
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = theta_star(c, k, X[i], &DF[i]);
/*
		printf("%02d:\ttheta*(%f) = %f\ttheta*'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Theta*(x)
	printf("Plotting theta*(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &F[1]);
	//	Show Theta*'(x)
	printf("Plotting theta*'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	// Test Theta - long range splitting function (finite)
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = theta(c, k, X[i], &DF[i]);
/*
		printf("%02d:\ttheta(%f) = %f\ttheta'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show Theta(x)
	printf("Plotting theta(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &F[1]);
	//	Show Theta'(x)
	printf("Plotting theta'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	// Test Gamma - long range splitting function (infinite)
	for (i = 1; i <= samples; i++)
	{
		X[i] = (2.0*(double)i/(double)samples);
		F[i] = gamma(c, k, X[i], &DF[i]);
/*
		printf("%02d:\tgamma(%f) = %f\tgamma'(%f) = %f\n",
				i, X[i], F[i], X[i], DF[i]);
*/
	}
	//	Show gamma(x)
	printf("Plotting gamma(x)  for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &F[1]);
	//	Show gamma'(x)
	printf("Plotting gamma'(x) for %2.1f <= x <= %2.1f...\t", X[1], X[samples]);
	plot2d(samples-1, &X[1], &DF[1]);

	/**************************************************************************/
	// Free allocated memory
	dynfree(c);
	dynfree(X);
	dynfree(F);
	dynfree(DF);
}

void stencil_naive(short p, double a, double h, short degree, double* omegap, short k, double* c, STENCIL* Ki, STENCIL* Kt)
{
	double		alpha = a/h;
	short		radius = (short)ceil(2*alpha);
	short		size = 2*radius+1;
	short		Ksize = size;//+2*degree;
	double*		Gamma = NULL;
	double*		KZ = NULL;
	double*		KY = NULL;
	double*		KX = NULL;
	long		x = 0;
	long		y = 0;
	long		z = 0;
	long		m = 0;
	double		d = 0.0;
	FILE*		fp = NULL;

	Gamma = (double*) dynvec(size*size*size, sizeof(double));
	KZ = (double*) dynvec(Ksize*Ksize*Ksize, sizeof(double));
	KY = (double*) dynvec(Ksize*Ksize*Ksize, sizeof(double));
	KX = (double*) dynvec(Ksize*Ksize*Ksize, sizeof(double));

	fp = fopen("Gamma.dat", "w");
	for (z = 0; z < size; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				d = h*sqrt((double)(x-radius)*(x-radius) + (double)(y-radius)*(y-radius) + (double)(z-radius)*(z-radius))/a;
				Gamma[z*size*size + y*size + x] = theta(c,k,d,NULL);
				fprintf(fp, "%+e ", Gamma[z*size*size + y*size + x]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KZ.dat", "w");
	//	Apply in Z direction
	for (z = 0; z < size; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				KZ[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[0]*Gamma[z*size*size + y*size + x];
				for (m = 1; m <= MIN(degree,size-1-z); m++)
				{
					KZ[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*Gamma[(z+m)*size*size + y*size + x];
				}
				for (m = 1; m <= MIN(degree,z); m++)
				{
					KZ[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*Gamma[(z-m)*size*size + y*size + x];
				}
				fprintf(fp, "%+e ", KZ[(z)*Ksize*Ksize + (y)*Ksize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KY.dat", "w");
	//	Apply in Y direction
	for (z = 0; z < size; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				KY[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[0]*KZ[z*size*size + y*size + x];
				for (m = 1; m <= MIN(degree,size-1-y,); m++)
				{
					KY[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*KZ[(z)*size*size + (y+m)*size + x];
				}
				for (m = 1; m <= MIN(degree,y); m++)
				{
					KY[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*KZ[(z)*size*size + (y-m)*size + x];
				}
				fprintf(fp, "%+e ", KY[(z)*Ksize*Ksize + (y)*Ksize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KX.dat", "w");
	//	Apply in X direction
	for (z = 0; z < size; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				KX[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[0]*KY[z*size*size + y*size + x];
				for (m = 1; m <= MIN(degree,size-1-x,); m++)
				{
					KX[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*KY[(z)*size*size + (y)*size + (x+m)];
				}
				for (m = 1; m <= MIN(degree,x); m++)
				{
					KX[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*KY[(z)*size*size + (y)*size + (x-m)];
				}
				fprintf(fp, "%+e ", KX[(z)*Ksize*Ksize + (y)*Ksize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	dynfree(KZ);
	dynfree(KY);
	dynfree(KX);
	dynfree(Gamma);
}

/*
Builds Gamma, anti-blurring operator, etc
*/
void test_preprocessing(void)
{
	short			k = 0;
	double			a = 0;
	double			h = 0.0;
	double			alpha = 0.0;
	short			p = 0;
	short			p_2 = 0;
	short			mu = 0;
	double*			c = NULL;
	double*			omegap = NULL;
	STENCIL			theta;
	STENCIL			gamma;
	STENCIL			g2g;
	STENCIL			tg2g;

	//	Get parameters from user
	printf("Please enter p: ");
	scanf("%hd", &p);

	printf("Please enter mu: ");
	scanf("%hd", &mu);

	printf("Enter k: ");
	scanf("%hd", &k);

	printf("Enter a: "); // alpha = a/h, radius = 2*alpha
	scanf("%lf", &a);

	printf("Enter h: "); // alpha = a/h, radius = 2*alpha
	scanf("%lf", &h);

	//	Compute computed values
	p_2 = p/2;
	alpha = a/h;

	// Dynamically allocate memory
	c = (double*) dynvec(k+1,sizeof(double));
	omegap = (double*) dynvec(p_2+mu+1, sizeof(double));

	//	Pre-pre-processing
	gamma_init(k, c);
	compute_omega_prime(p, mu, omegap);

	//	Pre-processing (Intermediate levels)
	stencil_initialize(&theta, (long) ceil(2.0*alpha), STENCIL_SHAPE_SPHERE);
	stencil_populate(&theta, c, k, STENCIL_FUNCTION_TYPE_THETA, h/a);
	stencil_display(&theta, h/a);

	stencil_initialize(&g2g, theta.size, theta.shape);
	stencil_shift(&theta, p_2 + mu, omegap, &g2g);
	stencil_display(&g2g, 1.0);

	//	Pre-processing (Top level)
	stencil_initialize(&gamma, 10 * (long) ceil(2.0*alpha), STENCIL_SHAPE_CUBE);
	stencil_populate(&gamma, c, k, STENCIL_FUNCTION_TYPE_GAMMA, h/a);
	//stencil_display(&gamma, h/a);

	stencil_initialize(&tg2g, gamma.size, gamma.shape);
	stencil_shift(&gamma, p_2 + mu, omegap, &tg2g);
	//stencil_display(&tg2g, 0.0);

	stencil_naive(p, a, h, p_2+mu, omegap, k, c, &g2g, &tg2g);

	//	Free dynamically allocated memory
	dynfree(c);
	dynfree(omegap);
	stencil_free(&theta);
	stencil_free(&gamma);
	stencil_free(&g2g);
	stencil_free(&tg2g);
}

// End of file