//-------|---------|---------|---------|---------|---------|---------|---------|
// all.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>

// TODO: reference additional headers your program requires here
#include <stdlib.h>
#include <assert.h>
#include <math.h>	//	fabs, ceil
#include <string.h> //	for memset and string functions
#include <mkl.h>

// Project specific
#include "memory.h"
#include "b_spline.h"
#include "even_powers.h"
#include "output.h"
#include "interpolant.h"

// Global constant(s):
#define MEM_ALIGN	64

// End of file