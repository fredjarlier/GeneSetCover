/**
 * \file mpi_globals.h
 * \author Frederic JARLIER
 */


#ifndef _MPI_GLOBALS_H_
#define _MPI_GLOBALS_H_

#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <errno.h>
#include <mpi.h>

extern int g_rank; // individual processor rank in MPI_COMM_WORLD
extern int g_size; // the total number of processors in MPI_COMM_WORLD
extern int g_iter; // what iterator are we on (used for writing files by each process)
extern FILE * g_fd; // global file descriptor used to write partial SAM file by each process
extern char g_name[MPI_MAX_PROCESSOR_NAME]; // name of node used, ie 'node9'
extern double g_start;
extern double g_end;
extern double g_elapsed;

#define DEBUG(format, ...) fprintf(stderr, "DEBUG %s:%d" format "\n", __FILE__, __LINE__, ##__VA_ARGS__);
#endif

/*
 	char tempString[128];
	char timestamp[100];
	struct timeval tv;
	struct tm *timeptr;
	gettimeofday(&tv, NULL);
	timeptr = localtime(&tv.tv_sec);

	strftime(timestamp, sizeof(timestamp), "%H:%M:%S", timeptr);

	int tempLength = snprintf(&tempString[0], 64, "[%s:%03d - P%d (%s) - %s] ", timestamp,

			(int)(tv.tv_usec / 1000), g_rank, g_name, __FUNCTION__);

	tempLength += snprintf(&tempString[tempLength], 64, format,__VA_ARGS__);
*/

extern void mpi_print_free_memory();
extern void mpi_open_psam_file();
extern void mpi_close_psam_file();
