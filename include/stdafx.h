//-------|---------|---------|---------|---------|---------|---------|---------|
// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

//#include "targetver.h" // WINDOWS SPECIFIC

#include <stdio.h>
//#include <tchar.h>	// WINDOWS SPECIFIC

// TODO: reference additional headers your program requires here
#include <stdlib.h>
#include <assert.h>
#include <math.h>	// fabs, ceil
#include <string.h> // for memset and string functions
#include <mkl.h>

// Project specific
#include "utility.h"
#include "phi.h"
#include "gamma.h"
#include "output.h"

// Global constant(s):
#define MEM_ALIGN	64

// End of file