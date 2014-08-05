//-------|---------|---------|---------|---------|---------|---------|---------|
/*
program output functions

Files needed:
	command_file	(may have generic template stored, or create in code)
	data_file		(must be passed into routines)
	plot_file		(gnuplot will create, but we will include if we create
					command file in code)
*/

#include "stdafx.h"

/*
Display 2D data only
	(uses auxiliary tmp file w/ standard I/O routines)
*/
void plot2d(int data_points, double* x, double* f)
{
	// Temporary data file name
	char* tmp_data_file = NULL;
	char* tmp_cmd_file = NULL;

	// Make sure incoming pointers are valid
	assert(data_points > 0);
	assert(x != NULL);
	assert(f != NULL);

	// Allocate memory for file name(s)
	tmp_data_file = (char*) dynvec(GP_DATA_DIR_LEN+GP_DATA_TMP_LEN+1,
									sizeof(char));
	tmp_cmd_file = (char*) dynvec(GP_CMD_DIR_LEN+GP_CMD_TEMPLATE_LEN+1,
									sizeof(char));

	// Set file name(s)
	sprintf(tmp_data_file, "%s%s", GP_DATA_DIR, GP_DATA_TMP);
	sprintf(tmp_cmd_file, "%s%s", GP_CMD_DIR, GP_CMD_TEMPLATE);

	// Save data to file
	save2d(data_points, x, f, tmp_data_file);

	// Call plotf
	plotf2d(tmp_cmd_file, tmp_data_file);

	// delete temporary data file
	if (remove(tmp_data_file))
	{
		printf("Error removing temporary file <%s> in plot\n", tmp_data_file);
	}

	// Free dynamically allocated memory
	dynfree(tmp_data_file);
	dynfree(tmp_cmd_file);
}

/*
Save 2D data and display 2D plot
*/
void plots2d(int data_points, double* x, double* f, char* save_file)
{
	// Command file name
	char* tmp_cmd_file = NULL;

	// Make sure incoming pointers are valid
	assert(save_file != NULL);
	assert(data_points > 0);
	assert(x != NULL);
	assert(f != NULL);

	// Allocate memory for file name(s)
	tmp_cmd_file = (char*) dynvec(GP_CMD_DIR_LEN+GP_CMD_TEMPLATE_LEN+1,
									sizeof(char));

	// Set file name(s)
	sprintf(tmp_cmd_file, "%s%s", GP_CMD_DIR, GP_CMD_TEMPLATE);

	// Save data to file
	save2d(data_points, x, f, save_file);

	// Call plotf
	plotf2d(tmp_cmd_file, save_file);

	// Free dynamically allocated memory
	dynfree(tmp_cmd_file);
}

/*
Save 2D data to file (cross-platform file I/O)
*/
void save2d(int data_points, double* x, double* f, char* save_file)
{
	FILE* fp = NULL;
	char buf[128] = {0};
	size_t bytes = 0;
	size_t buflen = 0;
	size_t sizeofchar = sizeof(char);
	int i = 0;

	// Make sure incoming pointers are valid
	assert(save_file != NULL);
	assert(data_points > 0);
	assert(x != NULL);
	assert(f != NULL);

	// Create temporary data file
	fp = fopen(save_file, "w");
	assert(fp != NULL);

	// Write contents of file
	for (i = 0; i < data_points; i++)
	{
		// FIXME - There is probably a better way to format this file... binary data would be most accurate, right?
		buflen = sprintf(buf, "%.32f %.32f\n", x[i], f[i]);
		bytes = fwrite(buf, sizeofchar, buflen, fp);
		if (bytes < buflen)
		{
			printf("<%d> bytes written to temporary file <%s> instead of <%d>\n", bytes, save_file, buflen);
			break;
		}
	}

	// Close file and flush IO stream to disk
	if (fclose(fp))
	{
		printf("Error closing temporary file in plot\n");
		return;
	}
}

/*
Driver for gnuplot to display data from data_file according to 
	format specified in command_file.
	(uses standard string manipulation functions)
	(uses standard system() process control function)
*/
void plotf2d(char* command_file, char* data_file)
{
	const char* cmdpre = "gnuplot ";
	const char* cmddata = "-e \"data_file='%s'\" ";
	char* cmd = NULL;
	int cmdlen = 0;
	int rc = 0;

	// Make sure that incoming pointers are valid
	assert(command_file != NULL);
	assert(data_file != NULL);

	// Compute length of string needed
	cmdlen += strlen(cmdpre);
	cmdlen += (strlen(cmddata)  - 2); // -2 for %s
	cmdlen += strlen(data_file);
	cmdlen += strlen(command_file);

	// Allocate memory for command string (+1 for NULL-terminator)
	cmd = (char*) dynvec(cmdlen+1, sizeof(char));

	// Copy the command prefix into the command string
	strncpy(cmd, cmdpre, strlen(cmdpre));
	// Copy the data file argument into the command string
	sprintf(&cmd[strlen(cmd)], cmddata, data_file);
	// Concatenate the filename to the command string
	strncat(cmd, command_file, cmdlen);

	// Have the system execute the command string
	if ((rc = system(cmd)))
	{
		// Display the return code from the command if error
		printf("gnuplot returned <%d>\n", rc);
	}

	// Release the allocated memory for the command string
	dynfree(cmd);
}

// End of file