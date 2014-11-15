//-------|---------|---------|---------|---------|---------|---------|---------|
/*
output.h - program output functions
*/

#ifndef	OUTPUT_H
#define OUTPUT_H

#include "all.h"
#include "memory.h"

#define GP_DATA_TMP			"tmp_data.dat"

/*
The follwoing preprocessor logic came from:
http://stackoverflow.com/questions/5919996/
how-to-detect-reliably-mac-os-x-ios-linux-windows-in-c-preprocessor
*/
#ifdef _WIN32
	//define something for Windows (32-bit and 64-bit, this part is common)
	#define GP_DATA_DIR			"../../../data/"
	#define GP_CMD_DIR			"../../../gnuplot/"
	#define GP_CMD_TEMPLATE		"tmp_data_template_windows.gp"
	#define GP_TERM				"windows"
	#ifdef _WIN64
		//define something for Windows (64-bit only)
	#endif
#elif __APPLE__
	#include "TargetConditionals.h"
	#if TARGET_IPHONE_SIMULATOR
		// iOS Simulator
	#elif TARGET_OS_IPHONE
		// iOS device
	#elif TARGET_OS_MAC
		// Other kinds of Mac OS
		#define GP_DATA_DIR			"data/"
		#define GP_CMD_DIR			"gnuplot/"
		#define GP_CMD_TEMPLATE		"tmp_data_template_aqua.gp"
		#define GP_TERM				"aqua"
	#else
		// Unsupported platform
	#endif
#elif __linux
	// linux
#elif __unix // all unices not caught above
	// Unix
#elif __posix
	// POSIX
#endif

#define	GP_DATA_TMP_LEN			strlen(GP_DATA_TMP)
#define	GP_CMD_TEMPLATE_LEN		strlen(GP_CMD_TEMPLATE)
#define	GP_DATA_DIR_LEN			strlen(GP_DATA_DIR)
#define	GP_CMD_DIR_LEN			strlen(GP_CMD_DIR)
#define	GP_TERM_LEN				strlen(GP_TERM)

#define GAMMA_DATA			"gamma.dat"
#define THETA_DATA			"theta.dat"
#define PHI_DATA			"phi.dat"
#define PHI_DATA_C1			"phiC1.dat"

#define GAMMA_DATA_LEN		strlen(GAMMA_DATA)
#define THETA_DATA_LEN		strlen(THETA_DATA)
#define PHI_DATA_LEN		strlen(PHI_DATA)
#define PHI_DATA_C1_LEN		strlen(PHI_DATA_C1)

void plot2d(int data_points, double* x, double* f);
void plots2d(int data_points, double* x, double* f, char* save_file);
void save2d(int data_points, double* x, double* f, char* save_file);
void plotf2d(char* command_file, char* data_file);

#endif

// End of file