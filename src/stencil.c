//-------|---------|---------|---------|---------|---------|---------|---------|
/*
stencil.c - 
*/

#include "stencil.h"

void stencil_initialize(STENCIL* s, long Size, short Shape)
{
	double			r = 0.0;
	double			rr = 0.0;
	long			z = 0;
	long			y = 0;
	long			zz = 0;
	long			yyzz = 0;
	long			zi_2d = 0;

	assert(s != NULL);

	//	Initialize stencil
	s->Shape = Shape;
	s->Size = Size;
	s->Data = (double*) dynvec(STENCIL_STORAGE(s->Size), sizeof(double));
	s->YMax = (long*) dynvec(s->Size+1, sizeof(long));						//	One max per z
	s->XMax = (long*) dynvec(STENCIL_STORAGE_2D(s->Size), sizeof(long));	//	One max per (y,z)

	//	Default all YMax and XMax values to -1
	for (z = 0; z <= s->Size; z++)
	{
		zi_2d = STENCIL_MAP_Y(z);
		s->YMax[z] = -1;
		for (y = 0; y <= s->YMax[z]; y++)
		{
			s->XMax[STENCIL_MAP_X(y) + zi_2d] = -1;
		}
	}

	//	Set up loop ranges
	if (s->Shape == STENCIL_SHAPE_SPHERE)
	{
		//	SPHERIC
		rr = (double) s->Size * s->Size;
		for (z = 0; z < s->Size; z++)	//	NOTE: Not going to s->Size
		{
			zz = z*z;
			zi_2d = STENCIL_MAP_Y(z);
			s->YMax[z] = MIN(z, (long) floor(sqrt(rr - zz)));
			for (y = 0; y <= s->YMax[z]; y++)
			{
				yyzz = y*y + zz;
				s->XMax[STENCIL_MAP_X(y) + zi_2d] = MIN(y, (long) floor(sqrt(rr - yyzz)));
			}
		}
	}
	else
	{
		//	CUBIC
		for (z = 0; z <= s->Size; z++)
		{
			zi_2d = STENCIL_MAP_Y(z);
			s->YMax[z] = z;
			for (y = 0; y <= s->YMax[z]; y++)
			{
				s->XMax[STENCIL_MAP_X(y) + zi_2d] = y;
			}
		}
	}
}

void stencil_populate(STENCIL* s, SOFTENER* Softener, short FunctionType, double Scale)
{
	void			(*f)(void*,long,double*,double*,double*);
	long			z = 0;
	long			y = 0;
	long			x = 0;
	long			zz = 0;
	long			zi = 0;
	long			zi_2d = 0;
	long			yyzz = 0;
	long			yi = 0;
	double			X = 0;
	double			DF = 0;

	assert(s != NULL);
	assert(Softener != NULL);

	if (FunctionType == STENCIL_FUNCTION_TYPE_THETA)
		f = Softener->split;
	else
		f = Softener->soften;

	for (z = 0; z <= s->Size; z++)
	{
		zz = z*z;
		zi = STENCIL_MAP_Z(z);
		zi_2d = STENCIL_MAP_Y(z);
		for (y = 0; y <= s->YMax[z]; y++)
		{
			yyzz = y*y + zz;
			yi = zi + STENCIL_MAP_Y(y);
			for (x = 0; x <= s->XMax[STENCIL_MAP_X(y) + zi_2d]; x++)
			{
				X = Scale*sqrt((double)x*x + yyzz);
				(*f)((void*)Softener, 1, &X, &s->Data[yi + STENCIL_MAP_X(x)], &DF);
			}
		}
	}
}

void stencil_display(STENCIL* s, double h_a)
{
	long				z = 0;
	long				y = 0;
	long				x = 0;
	long				idx = 0;
	double				d = 0.0;

	for (z = 0; z <= s->Size; z++)
	{
		for (y = 0; y <= s->YMax[z]; y++)
		{
			for (x = 0; x <= s->XMax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
				d = h_a*sqrt((double) x*x + y*y + z*z);
				idx = STENCIL_MAP_X(x) + STENCIL_MAP_Y(y) + STENCIL_MAP_Z(z);
				printf("(%ld,%ld,%ld) -> %e -> %+e\n", x,y,z, d, s->Data[idx]);
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
	long				r = s->Size;
	double*				tmp1 = NULL;
	double*				tmp2 = NULL;

	assert(s != NULL);
	assert(K != NULL);
	assert(s->Size == K->Size);

	//	Create memory for intermediate stencils
	tmp1 = (double*) dynvec((K->Size+1)*(K->Size+1)*(K->Size+2)/2,sizeof(double));

	//	Apply anti-blurring operator to s in Z direction, i.e., (A_z)s
//***NOTE: COULD RESTRICT LOOPS TO SPHERIC INDEXES OF K, I THINK***
	for (z = 0; z <= K->Size; z++)
	{
		zz = STENCIL_MAP_Z2(K->Size,z);
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= y; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->Data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				for (m = 1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-y,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z-m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = z-y+1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = z+x+1; m <= MIN(y+z,n); m++)
				{
                    tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}
		}

		for (y = z+1; y <= K->Size; y++)
		{
			yy = STENCIL_MAP_Y2(y);
			for (x = 0; x <= z; x++)
			{
// x <= z < y
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z)+STENCIL_MAP_X(x)];
				for (m = 1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}

			for (x = z+1; x <= y; x++)
			{
// z < x <= y
				xx = STENCIL_MAP_X2(x);
				tmp1[zz+yy+xx] = omegap[0]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z)];
				for (m = 1; m <= MIN(x-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z+m)];
				}
				for (m = x-z+1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= MIN(r-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= MIN(z+r,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}
		}
	}
/*
	//	Display "stacked" half-plane where x-y plane is symmetric
	printf("tmp1:\n");
	for (z = 0; z <= K->Size; z++)
	{
		zz = STENCIL_MAP_Z2(K->Size,z);//z*(K->Size+2)*(K->Size+1)/2;
		for (y = 0; y <= K->Size; y++)
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
	tmp2 = (double*) dynvec((K->Size+1)*(K->Size+1)*(K->Size+2)/2,sizeof(double));
	r = K->Size;	//	This represents the length of the "stacked" symmetric plane stencils
	for (z = 0; z <= K->Size; z++)
	{
		zz = STENCIL_MAP_Y2(z);
		for (y = 0; y <= z; y++)
		{
			yy = STENCIL_MAP_X2(y);
			for (x = 0; x <= y; x++)
			{
//	x <= y <= z
				xx = STENCIL_MAP_Z2(K->Size,x);
				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(y)+STENCIL_MAP_X2(x)];
				for (m = 1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y-x,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(y-m)+STENCIL_MAP_X2(x)];
				}
				for (m = y-x+1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}

			for (x = y+1; x <= z; x++)
			{
//	y <= x <= z
				xx = STENCIL_MAP_Z2(K->Size,x);
				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(x-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
				}
				for (m = x-y+1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}

			for (x = z+1; x <= K->Size; x++)
			{
//	y <= z <= x
				xx = STENCIL_MAP_Z2(K->Size,x);
				tmp2[xx+zz+yy] = omegap[0]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(x-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y+m)];
				}
				for (m = x-y+1; m <= MIN(r-y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(y+m)+STENCIL_MAP_X2(x)];
				}
				for (m = 1; m <= MIN(y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(y-m)];
				}
				for (m = y+1; m <= MIN(x+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(x)+STENCIL_MAP_X2(m-y)];
				}
				for (m = x+y+1; m <= MIN(r+y,n); m++)
				{
					tmp2[xx+zz+yy] += omegap[m]*tmp1[STENCIL_MAP_Z2(K->Size,z)+STENCIL_MAP_Y2(m-y)+STENCIL_MAP_X2(x)];
				}
			}
		}
	}

	//	Free dynamically allocated memory for tmp1
	dynfree(tmp1);
/*
	//	Display "stacked" half-plane where y-z plane is symmetric
	printf("\ntmp2\n");
	for (x = 0; x <= K->Size; x++)
	{
		xx = STENCIL_MAP_Z2(K->Size,x);
		for (z = 0; z <= K->Size; z++)
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
	r = K->Size;
	for (z = 0; z <= K->Size; z++)
	{
		zz = STENCIL_MAP_Z(z);
		for (y = 0; y <= K->YMax[z]; y++)
		{
			yy = STENCIL_MAP_Y(y);
			for (x = 0; x <= K->XMax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X(x);
				K->Data[zz+yy+xx] = omegap[0]*tmp2[STENCIL_MAP_Z2(K->Size,x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= MIN(r-x,n); m++)
				{
					K->Data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->Size,x+m)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = 1; m <= MIN(x,n); m++)
				{
					K->Data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->Size,x-m)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = x+1; m <= MIN(r+x,n); m++)
				{
					K->Data[zz+yy+xx] += omegap[m]*tmp2[STENCIL_MAP_Z2(K->Size,m-x)+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
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
	assert(s->Size - K->Size == degree);

	lradius = s->Size;	//	large radius
	sradius = K->Size;	//	small radius
/*
Gamma:	[l,l,l]	-> s
KZ:		[l,l,s] -> tmp1 (store as [l,l,s] = [x,y,z])
KY:		[l,s,s] -> tmp2 (store as [s,s,l] = [y,z,x])
KX:		[s,s,s]	-> K
NOTE: l[arge] = s->Size >= K->Size = s[mall]
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
				tmp1[zz+yy+xx] = omegap[0]*s->Data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				for (m = 1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-y,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z-m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = z-y+1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = z+x+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
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
				tmp1[zz+yy+xx] = omegap[0]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z)+STENCIL_MAP_X(x)];
				for (m = 1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z-x,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z-m)+STENCIL_MAP_X(x)];
				}
				for (m = z-x+1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}

			for (x = z+1; x <= y; x++)
			{
				xx = STENCIL_MAP_X2(x);
// z < x <= y
				tmp1[zz+yy+xx] = omegap[0]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z)];
				for (m = 1; m <= MIN(x-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z+m)];
				}
				for (m = x-z+1; m <= MIN(y-z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(z+m)+STENCIL_MAP_X(x)];
				}
				for (m = y-z+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(z+m)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
				for (m = 1; m <= MIN(z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(z-m)];
				}
				for (m = z+1; m <= MIN(x+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(x)+STENCIL_MAP_X(m-z)];
				}
				for (m = x+z+1; m <= MIN(y+z,n); m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(y)+STENCIL_MAP_Y(m-z)+STENCIL_MAP_X(x)];
				}
				for (m = z+y+1; m <= n; m++)
				{
					tmp1[zz+yy+xx] += omegap[m]*s->Data[STENCIL_MAP_Z(m-z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)];
				}
			}
		}
	}
/*
	//	Display "stacked" half-plane where x-y plane is symmetric
	printf("tmp1:\n");
	for (z = 0; z <= K->Size; z++)
	{
		zz = z*(lradius+1)*(lradius+2)/2;
		for (y = 0; y <= K->Size; y++)
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
	for (x = 0; x <= K->Size; x++)
	{
		xx = x*(sradius+1)*(sradius+2)/2;
		for (z = 0; z <= K->Size; z++)
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
		for (y = 0; y <= K->YMax[z]; y++)
		{
			yy = STENCIL_MAP_Y(y);
			for (x = 0; x <= K->XMax[STENCIL_MAP_X(y) + STENCIL_MAP_Y(z)]; x++)
			{
// x <= y <= z
				xx = STENCIL_MAP_X(x);
				K->Data[zz+yy+xx] = omegap[0]*tmp2[(x)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				for (m = 1; m <= n; m++)
				{
					K->Data[zz+yy+xx] += omegap[m]*tmp2[(x+m)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = 1; m <= MIN(x,n); m++)
				{
					K->Data[zz+yy+xx] += omegap[m]*tmp2[(x-m)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
				for (m = x+1; m <= n; m++)
				{
					K->Data[zz+yy+xx] += omegap[m]*tmp2[(m-x)*(sradius+1)*(sradius+2)/2+STENCIL_MAP_Y2(z)+STENCIL_MAP_X2(y)];
				}
			}
		}
	}

	//	Free dynamically allocated memory for tmp2
	dynfree(tmp2);
}

void stencil_free(STENCIL* s)
{
	if (s->Data != NULL)
		dynfree(s->Data);

	if (s->XMax != NULL)
		dynfree(s->XMax);

	if (s->YMax != NULL)
		dynfree(s->YMax);
}

#if 0
void stencil_naive(short p, double a, double h, short degree, double* omegap, short k, double* c, STENCIL* Ki)
{
	double		alpha = a/h;
	short		radius = (short)ceil(2*alpha);
	short		Size = 2*radius+1;
	short		KSize = Size;//+2*degree;
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

	Gamma = (double*) dynvec(Size*Size*Size, sizeof(double));
	KZ = (double*) dynvec(KSize*KSize*KSize, sizeof(double));
	KY = (double*) dynvec(KSize*KSize*KSize, sizeof(double));
	KX = (double*) dynvec(KSize*KSize*KSize, sizeof(double));

	fp = fopen("Gamma.dat", "w");
	for (z = 0; z < Size; z++)
	{
		for (y = 0; y < Size; y++)
		{
			for (x = 0; x < Size; x++)
			{
				d = h*sqrt((double)(x-radius)*(x-radius) + (double)(y-radius)*(y-radius) + (double)(z-radius)*(z-radius))/a;
				Gamma[z*Size*Size + y*Size + x] = theta(c,k,d,NULL);
				fprintf(fp, "%+e ", Gamma[z*Size*Size + y*Size + x]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KZ.dat", "w");
	//	Apply in Z direction
	for (z = 0; z < Size; z++)
	{
		for (y = 0; y < Size; y++)
		{
			for (x = 0; x < Size; x++)
			{
				KZ[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[0]*Gamma[z*Size*Size + y*Size + x];
				for (m = 1; m <= MIN(degree,Size-1-z); m++)
				{
					KZ[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*Gamma[(z+m)*Size*Size + y*Size + x];
				}
				for (m = 1; m <= MIN(degree,z); m++)
				{
					KZ[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*Gamma[(z-m)*Size*Size + y*Size + x];
				}
				fprintf(fp, "%+e ", KZ[(z)*KSize*KSize + (y)*KSize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KY.dat", "w");
	//	Apply in Y direction
	for (z = 0; z < Size; z++)
	{
		for (y = 0; y < Size; y++)
		{
			for (x = 0; x < Size; x++)
			{
				KY[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[0]*KZ[z*Size*Size + y*Size + x];
				for (m = 1; m <= MIN(degree,Size-1-y); m++)
				{
					KY[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*KZ[(z)*Size*Size + (y+m)*Size + x];
				}
				for (m = 1; m <= MIN(degree,y); m++)
				{
					KY[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*KZ[(z)*Size*Size + (y-m)*Size + x];
				}
				fprintf(fp, "%+e ", KY[(z)*KSize*KSize + (y)*KSize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KX.dat", "w");
	//	Apply in X direction
	for (z = 0; z < Size; z++)
	{
		for (y = 0; y < Size; y++)
		{
			for (x = 0; x < Size; x++)
			{
				KX[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[0]*KY[z*Size*Size + y*Size + x];
				for (m = 1; m <= MIN(degree,Size-1-x); m++)
				{
					KX[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*KY[(z)*Size*Size + (y)*Size + (x+m)];
				}
				for (m = 1; m <= MIN(degree,x); m++)
				{
					KX[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*KY[(z)*Size*Size + (y)*Size + (x-m)];
				}
				fprintf(fp, "%+e ", KX[(z)*KSize*KSize + (y)*KSize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (z = 0; z <= Ki->Size; z++)
	{
		for (y = 0; y <= Ki->YMax[z]; y++)
		{
			for (x = 0; x <= Ki->XMax[STENCIL_MAP_Y(z)+STENCIL_MAP_X(y)]; x++)
			{
				// +radius for KX
				err = fabs(Ki->Data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)] - KX[(z+radius)*KSize*KSize+(y+radius)*KSize+(x+radius)]) / fabs(KX[(z+radius)*KSize*KSize+(y+radius)*KSize+(x+radius)]);
//				printf("(%02ld,%02ld,%02ld) Ki = %+e, KX = %+e, err = %+e\n", x,y,z, Ki->Data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)], KX[(z+radius)*KSize*KSize+(y+radius)*KSize+(x+radius)], err);
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
	short		radius = Kt->Size+degree;
	short		Size = 2*Kt->Size+2*degree+1;
	short		KSize = 2*Kt->Size+1;
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
	Gamma = (double*) dynvec(Size*Size*Size, sizeof(double));
	KZ = (double*) dynvec(Size*Size*KSize, sizeof(double));
	KY = (double*) dynvec(Size*KSize*KSize, sizeof(double));
	KX = (double*) dynvec(KSize*KSize*KSize, sizeof(double));

	fp = fopen("Gamma.dat", "w");
	for (z = 0; z < Size; z++)
	{
		for (y = 0; y < Size; y++)
		{
			for (x = 0; x < Size; x++)
			{
				d = h*sqrt((double)(x-radius)*(x-radius) + (double)(y-radius)*(y-radius) + (double)(z-radius)*(z-radius))/a;
				Gamma[z*Size*Size + y*Size + x] = _gamma(c,k,d,NULL);
				fprintf(fp, "%+e ", Gamma[z*Size*Size + y*Size + x]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KZ.dat", "w");
	//	Apply in Z direction
	for (z = 0; z < KSize; z++)
	{
		for (y = 0; y < Size; y++)
		{
			for (x = 0; x < Size; x++)
			{
				KZ[(z)*Size*Size + (y)*Size + (x)] += omegap[0]*Gamma[(z+degree)*Size*Size + (y)*Size + (x)];
				for (m = 1; m <= degree; m++)
				{
					KZ[(z)*Size*Size + (y)*Size + (x)] += omegap[m]*Gamma[(z+degree+m)*Size*Size + (y)*Size + (x)];
				}
				for (m = 1; m <= degree; m++)
				{
					KZ[(z)*Size*Size + (y)*Size + (x)] += omegap[m]*Gamma[(z+degree-m)*Size*Size + (y)*Size + (x)];
				}
				fprintf(fp, "%+e ", KZ[(z)*Size*Size + (y)*Size + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KY.dat", "w");
	//	Apply in Y direction
	for (z = 0; z < KSize; z++)
	{
		for (y = 0; y < KSize; y++)
		{
			for (x = 0; x < Size; x++)
			{
				KY[(z)*KSize*Size + (y)*Size + (x)] += omegap[0]*KZ[z*Size*Size + (y+degree)*Size + x];
				for (m = 1; m <= degree; m++)
				{
					KY[(z)*KSize*Size + (y)*Size + (x)] += omegap[m]*KZ[(z)*Size*Size + (y+degree+m)*Size + x];
				}
				for (m = 1; m <= degree; m++)
				{
					KY[(z)*KSize*Size + (y)*Size + (x)] += omegap[m]*KZ[(z)*Size*Size + (y+degree-m)*Size + x];
				}
				fprintf(fp, "%+e ", KY[(z)*KSize*Size + (y)*Size + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("KX.dat", "w");
	//	Apply in X direction
	for (z = 0; z < KSize; z++)
	{
		for (y = 0; y < KSize; y++)
		{
			for (x = 0; x < KSize; x++)
			{
				KX[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[0]*KY[z*KSize*Size + y*Size + (x+degree)];
				for (m = 1; m <= degree; m++)
				{
					KX[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*KY[(z)*KSize*Size + (y)*Size + (x+degree+m)];
				}
				for (m = 1; m <= degree; m++)
				{
					KX[(z)*KSize*KSize + (y)*KSize + (x)] += omegap[m]*KY[(z)*KSize*Size + (y)*Size + (x+degree-m)];
				}
				fprintf(fp, "%+e ", KX[(z)*KSize*KSize + (y)*KSize + (x)]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (z = 0; z <= Kt->Size; z++)
	{
		for (y = 0; y <= Kt->YMax[z]; y++)
		{
			for (x = 0; x <= Kt->XMax[STENCIL_MAP_Y(z)+STENCIL_MAP_X(y)]; x++)
			{
				// +radius for KX
				err = fabs(Kt->Data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)] - KX[(z+Kt->Size)*KSize*KSize+(y+Kt->Size)*KSize+(x+Kt->Size)]) / fabs(KX[(z+Kt->Size)*KSize*KSize+(y+Kt->Size)*KSize+(x+Kt->Size)]);
//				printf("(%02ld,%02ld,%02ld) Kt = %+e, KX = %+e, err = %+e\n", x,y,z, Kt->Data[STENCIL_MAP_Z(z)+STENCIL_MAP_Y(y)+STENCIL_MAP_X(x)], KX[(z+Kt->Size)*KSize*KSize+(y+Kt->Size)*KSize+(x+Kt->Size)], err);
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
#endif

//	End of file