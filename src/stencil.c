//-------|---------|---------|---------|---------|---------|---------|---------|
/*
stencil.c - 
*/

#include "stencil.h"

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
		f = &_gamma;

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
//***NOTE: COULD RESTRICT LOOPS TO SPHERIC INDEXES OF K, I THINK***
	for (z = 0; z <= K->size; z++)
	{
		zz = STENCIL_MAP_Z2(K->size,z);
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= y; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				for (m = 1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-y,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z-m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = z-y+1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = z+x+1; m <= MIN(y+z,n); m++)
				{
                    tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}
		}

		for (y = z+1; y <= K->size; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= z; x++)
			{
// x <= z < y
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z)+STENCIL_MAP_X(x)];
				for (m = 1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}

			for (x = z+1; x <= y; x++)
			{
// z < x <= y
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z)];
				for (m = 1; m <= MIN(x-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z+m)];
				}
				for (m = x-z+1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
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
	//	Apply anti-blurring operator to (A_z)s in Y direction, i.e., (A_y)(A_z)s
	tmp2 = (double*) dynvec((K->size+1)*(K->size+1)*(K->size+2)/2,sizeof(double));
	r = K->size;	//	This represents the length of the "stacked" symmetric plane stencils
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
				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y)+STENCIL_MAP_X2(x)];
				for (m = 1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y-x,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y-m)+STENCIL_MAP_X2(x)];
				}
				for (m = y-x+1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}

			for (x = y+1; x <= z; x++)
			{
//	y <= x <= z
				xx = STENCIL_MAP_Z2(K->size,x);
				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(x-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
				}
				for (m = x-y+1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}

			for (x = z+1; x <= K->size; x++)
			{
//	y <= z <= x
				xx = STENCIL_MAP_Z2(K->size,x);
				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(x-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
				}
				for (m = x-y+1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
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
		xx = STENCIL_MAP_Z2(K->size,x);
		for (z = 0; z <= K->size; z++)
		{
			zz = STENCIL_MAP_Y2(z);
			for (y = 0; y <= z; y++)
			{
				yy = STENCIL_MAP_X2(y);
				printf("%+e ", tmp2[xx+zz+yy]);
			}
			printf("\n");
		}
		printf("\n");
	}
*/
	//	Apply anti-blurring operator to (A_y)(A_z)s in X direction, i.e., (A_x)(A_y)(A_z)s
	r = K->size;
	for (z = 0; z <= K->size; z++)
	{
		zz = STENCIL_MAP_Z(z);
		for (y = 0; y <= K->ymax[z]; y++)
		{
			yy = STENCIL_MAP_Y(y);
			for (x = 0; x <= K->xmax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X(x);
				K->data[zz+yy+xx] = omegap[0]*tmp2[STENCIL_MAP_Z2(K->size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(r-x,n); m++)
				{
					K->data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->size,x+m)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = 1; m <= MIN(x,n); m++)
				{
					K->data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->size,x-m)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = x+1; m <= MIN(r+x,n); m++)
				{
					K->data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->size,m-x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
			}
		}
	}

	//	Free dynamically allocated memory for tmp2
	dynfree(tmp2);
}

void stencil_shift_infinite(STENCIL* s, short degree, double* omegap, STENCIL* K)
{
	long				z = 0;
	long				y = 0;
	long				x = 0;
	short				n = degree;
	long				xx = 0, yy = 0, zz = 0;
	long				m = 0;
	double*				tmp1 = NULL;
	double*				tmp2 = NULL;
	long				lradius = 0;
	long				sradius = 0;

	assert(s != NULL);
	assert(K != NULL);
	assert(s->size - K->size == degree);

	lradius = s->size;	//	large radius
	sradius = K->size;	//	small radius
/*
Gamma:	[l,l,l]	-> s
KZ:		[l,l,s] -> tmp1 (store as [l,l,s] = [x,y,z])
KY:		[l,s,s] -> tmp2 (store as [s,s,l] = [y,z,x])
KX:		[s,s,s]	-> K
NOTE: l[arge] = s->size >= K->size = s[mall]
*/
	//	Create memory for intermediate stencils
	tmp1 = (double*) dynvec((sradius+1)*(lradius+1)*(lradius+2)/2,sizeof(double));

	//	Apply anti-blurring operator to s in Z direction, i.e., (A_z)s
//***NOTE: COULD RESTRICT LOOPS TO SPHERIC INDEXES OF K, I THINK***
	for (z = 0; z <= sradius; z++)
	{
		zz = z*(lradius+1)*(lradius+2)/2;
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= y; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				for (m = 1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-y,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z-m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = z-y+1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = z+x+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}
		}

		for (y = z+1; y <= lradius; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= z; x++)
			{
				xx = STENCIL_MAP_X2(x);
// x <= z < y
				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z)+STENCIL_MAP_X(x)];
				for (m = 1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}

			for (x = z+1; x <= y; x++)
			{
				xx = STENCIL_MAP_X2(x);
// z < x <= y
				tmp1[zz+yy+xx] = omegap[0]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z)];
				for (m = 1; m <= MIN(x-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z+m)];
				}
				for (m = x-z+1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}
		}
	}
/*
	//	Display "stacked" half-plane where x-y plane is symmetric
	printf("tmp1:\n");
	for (z = 0; z <= K->size; z++)
	{
		zz = z*(lradius+1)*(lradius+2)/2;
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
	//	Apply anti-blurring operator to (A_z)s in Y direction, i.e., (A_y)(A_z)s
	tmp2 = (double*) dynvec((lradius+1)*(sradius+1)*(sradius+2)/2,sizeof(double));
	for (z = 0; z <= sradius; z++)
	{
		zz = STENCIL_MAP_Y2(z);
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_X2(y);
			for (x = 0; x <= y; x++)
			{
//	x <= y <= z
				xx = x*(sradius+1)*(sradius+2)/2;
				tmp2[xx+zz+yy] = omegap[0]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(y)+STENCIL_MAP_X2(x)];
				for (m = 1; m <= n; m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y-x,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(y-m)+STENCIL_MAP_X2(x)];
				}
				for (m = y-x+1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= n; m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}

			for (x = y+1; x <= z; x++)
			{
//	y <= x <= z
				xx = x*(sradius+1)*(sradius+2)/2;
				tmp2[xx+zz+yy] = omegap[0]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(x-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
				}
				for (m = x-y+1; m <= n; m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= n; m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}

			for (x = z+1; x <= lradius; x++)
			{
//	y <= z <= x
				xx = x*(sradius+1)*(sradius+2)/2;
				tmp2[xx+zz+yy] = omegap[0]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(x-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
				}
				for (m = x-y+1; m <= n; m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= n; m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[z*(lradius+1)*(lradius+2)/2+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
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
		xx = x*(sradius+1)*(sradius+2)/2;
		for (z = 0; z <= K->size; z++)
		{
			zz = STENCIL_MAP_Y2(z);
			for (y = 0; y <= z; y++)
			{
				yy = STENCIL_MAP_X2(y);
				printf("%+e ", tmp2[xx+zz+yy]);
			}
			printf("\n");
		}
		printf("\n");
	}
*/
	//	Apply anti-blurring operator to (A_y)(A_z)s in X direction, i.e., (A_x)(A_y)(A_z)s
	for (z = 0; z <= sradius; z++)
	{
		zz = STENCIL_MAP_Z(z);
		for (y = 0; y <= K->ymax[z]; y++)
		{
			yy = STENCIL_MAP_Y(y);
			for (x = 0; x <= K->xmax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X(x);
				K->data[zz+yy+xx] = omegap[0]*tmp2[(x)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= n; m++)
				{
					K->data[zz+yy+xx] += omegap[m]*tmp2[(x+m)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = 1; m <= MIN(x,n); m++)
				{
					K->data[zz+yy+xx] += omegap[m]*tmp2[(x-m)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = x+1; m <= n; m++)
				{
					K->data[zz+yy+xx] += omegap[m]*tmp2[(m-x)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
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

void stencil_naive(short p, double a, double h, short degree, double* omegap, short k, double* c, STENCIL* Ki)
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
	double		minerr = 1.0;
	double		maxerr = 0.0;
	double		err = 0.0;

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
				for (m = 1; m <= MIN(degree,size-1-y); m++)
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
				for (m = 1; m <= MIN(degree,size-1-x); m++)
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

	for (z = 0; z <= Ki->size; z++)
	{
		for (y = 0; y <= Ki->ymax[z]; y++)
		{
			for (x = 0; x <= Ki->xmax[STENCIL_MAP_Y(z)+STENCIL_MAP_X(y)]; x++)
			{
				// +radius for KX
				err = fabs(Ki->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)] - KX[(z+radius)*Ksize*Ksize+(y+radius)*Ksize+(x+radius)]) / fabs(KX[(z+radius)*Ksize*Ksize+(y+radius)*Ksize+(x+radius)]);
//				printf("(%02ld,%02ld,%02ld) Ki = %+e, KX = %+e, err = %+e\n", x,y,z, Ki->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)], KX[(z+radius)*Ksize*Ksize+(y+radius)*Ksize+(x+radius)], err);
				if (err > maxerr)
					maxerr = err;
				if (err < minerr)
					minerr = err;
			}
		}
	}
	printf("Ki: minerr = %+e, maxerr = %+e\n", minerr, maxerr);

	dynfree(KZ);
	dynfree(KY);
	dynfree(KX);
	dynfree(Gamma);
}

void stencil_naive_top(short p, double a, double h, short degree, double* omegap, short k, double* c, STENCIL* Kt)
{
	short		radius = Kt->size+degree;
	short		size = 2*Kt->size+2*degree+1;
	short		Ksize = 2*Kt->size+1;
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
	double		minerr = 1.0;
	double		maxerr = 0.0;
	double		err = 0.0;
/*
Gamma:	[l, l, l]
KZ:		[l, l, s]
KY:		[l, s, s]
KX:		[s, s, s]
*/
	Gamma = (double*) dynvec(size*size*size, sizeof(double));
	KZ = (double*) dynvec(size*size*Ksize, sizeof(double));
	KY = (double*) dynvec(size*Ksize*Ksize, sizeof(double));
	KX = (double*) dynvec(Ksize*Ksize*Ksize, sizeof(double));

	fp = fopen("Gamma.dat", "w");
	for (z = 0; z < size; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				d = h*sqrt((double)(x-radius)*(x-radius) + (double)(y-radius)*(y-radius) + (double)(z-radius)*(z-radius))/a;
				Gamma[z*size*size + y*size + x] = _gamma(c,k,d,NULL);
				fprintf(fp, "%+e ", Gamma[z*size*size + y*size + x]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KZ.dat", "w");
	//	Apply in Z direction
	for (z = 0; z < Ksize; z++)
	{
		for (y = 0; y < size; y++)
		{
			for (x = 0; x < size; x++)
			{
				KZ[(z)*size*size + (y)*size + (x)] += omegap[0]*Gamma[(z+degree)*size*size + (y)*size + (x)];
				for (m = 1; m <= degree; m++)
				{
					KZ[(z)*size*size + (y)*size + (x)] += omegap[m]*Gamma[(z+degree+m)*size*size + (y)*size + (x)];
				}
				for (m = 1; m <= degree; m++)
				{
					KZ[(z)*size*size + (y)*size + (x)] += omegap[m]*Gamma[(z+degree-m)*size*size + (y)*size + (x)];
				}
				fprintf(fp, "%+e ", KZ[(z)*size*size + (y)*size + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KY.dat", "w");
	//	Apply in Y direction
	for (z = 0; z < Ksize; z++)
	{
		for (y = 0; y < Ksize; y++)
		{
			for (x = 0; x < size; x++)
			{
				KY[(z)*Ksize*size + (y)*size + (x)] += omegap[0]*KZ[z*size*size + (y+degree)*size + x];
				for (m = 1; m <= degree; m++)
				{
					KY[(z)*Ksize*size + (y)*size + (x)] += omegap[m]*KZ[(z)*size*size + (y+degree+m)*size + x];
				}
				for (m = 1; m <= degree; m++)
				{
					KY[(z)*Ksize*size + (y)*size + (x)] += omegap[m]*KZ[(z)*size*size + (y+degree-m)*size + x];
				}
				fprintf(fp, "%+e ", KY[(z)*Ksize*size + (y)*size + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KX.dat", "w");
	//	Apply in X direction
	for (z = 0; z < Ksize; z++)
	{
		for (y = 0; y < Ksize; y++)
		{
			for (x = 0; x < Ksize; x++)
			{
				KX[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[0]*KY[z*Ksize*size + y*size + (x+degree)];
				for (m = 1; m <= degree; m++)
				{
					KX[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*KY[(z)*Ksize*size + (y)*size + (x+degree+m)];
				}
				for (m = 1; m <= degree; m++)
				{
					KX[(z)*Ksize*Ksize + (y)*Ksize + (x)] += omegap[m]*KY[(z)*Ksize*size + (y)*size + (x+degree-m)];
				}
				fprintf(fp, "%+e ", KX[(z)*Ksize*Ksize + (y)*Ksize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (z = 0; z <= Kt->size; z++)
	{
		for (y = 0; y <= Kt->ymax[z]; y++)
		{
			for (x = 0; x <= Kt->xmax[STENCIL_MAP_Y(z)+STENCIL_MAP_X(y)]; x++)
			{
				// +radius for KX
				err = fabs(Kt->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)] - KX[(z+Kt->size)*Ksize*Ksize+(y+Kt->size)*Ksize+(x+Kt->size)]) / fabs(KX[(z+Kt->size)*Ksize*Ksize+(y+Kt->size)*Ksize+(x+Kt->size)]);
//				printf("(%02ld,%02ld,%02ld) Kt = %+e, KX = %+e, err = %+e\n", x,y,z, Kt->data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)], KX[(z+Kt->size)*Ksize*Ksize+(y+Kt->size)*Ksize+(x+Kt->size)], err);
				if (err > maxerr)
					maxerr = err;
				if (err < minerr)
					minerr = err;
			}
		}
	}
	printf("Kt: minerr = %+e, maxerr = %+e\n", minerr, maxerr);

	dynfree(KZ);
	dynfree(KY);
	dynfree(KX);
	dynfree(Gamma);
}

//	End of file