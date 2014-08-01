//-------|---------|---------|---------|---------|---------|---------|---------|
/*
program output functions
*/

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