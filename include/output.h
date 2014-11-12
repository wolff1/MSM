//-------|---------|---------|---------|---------|---------|---------|---------|
/*
output.h - program output functions
*/

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

/*
Display 2D data only
	(uses auxiliary tmp file w/ standard I/O routines)
*/
void plot2d(int data_points, double* x, double* f);

/*
Save 2D data and display 2D plot
*/
void plots2d(int data_points, double* x, double* f, char* save_file);

/*
Save 2D data to file (cross-platform file I/O)
*/
void save2d(int data_points, double* x, double* f, char* save_file);

/*
Driver for gnuplot to display data from data_file according to 
	format specified in command_file.
	(uses standard string manipulation functions)
	(uses standard system() process control function)
*/
void plotf2d(char* command_file, char* data_file);

// End of file