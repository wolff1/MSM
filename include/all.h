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

#endif

// End of file