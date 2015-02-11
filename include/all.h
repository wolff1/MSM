//-------|---------|---------|---------|---------|---------|---------|---------|
// all.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef ALL_H
#define ALL_H

#pragma once

#include <stdio.h>

// TODO: reference additional headers your program requires here
#include <stdlib.h>
#include <assert.h>
#include <math.h>	//	fabs, ceil
#include <string.h> //	for memset and string functions
#include <mkl.h>

// Global macros:
#define	MAX(x,y)						((x > y) ? (x) : (y))
#define	MIN(x,y)						((x < y) ? (x) : (y))
#define	PI								3.14159265358979323846
#define	IDX(x,y,z,XLength,YLength)		(z)*(XLength)*(YLength) + (y)*(XLength) + (x)	//	3D to 1D

#define	SORT_RTN_STNCL_IDX(x,y,z,q,m,n,idx) \
		if (x < y)	{if (x < z) {q = x;	if (y < z) {m = y; n = z;} else {m = z; n = y;}} else {q = z; m = x; n = y;}} \
		else		{if (y < z) {q = y;	if (x < z) {m = x; n = z;} else {m = z; n = x;}} else {q = z; m = y; n = x;}} \
		idx = STENCIL_MAP_X(q) + STENCIL_MAP_Y(m) + STENCIL_MAP_Z(n);

#endif

// End of file