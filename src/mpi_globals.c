/**
 * \file mpi_globals.c
 * \author Frederic JARLIER
 */

#include <errno.h>
#include <stdio.h>
#include <limits.h>
#include "mpi_globals.h"
#include "mpi.h"

int g_rank;
int g_size;
int g_iter = 0;

char g_name[MPI_MAX_PROCESSOR_NAME];

FILE * g_fd = NULL;

double g_start;
double g_end;
double g_elapsed = 0; // time elapised measurement

void mpi_print_free_memory(void)
{
	FILE * fp = popen("free -m | grep Mem | awk '{print $2}'", "r");
	if (NULL == fp)
	{
		// errno 12 means not enough memory
		DEBUG("COULD NOT READ FREE MEMORY, error=%d\n", errno);
		return;
	}
	char path[1035];
	while (fgets(path, sizeof(path)-1, fp) != NULL) {
		DEBUG(" \n Free MB: %s \n", path);
	}

	DEBUG("\n its done!! \n");
	pclose(fp);
}



int MPI_Large_bcast(void * buffer, unsigned long count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
	while (count > INT_MAX)
	{
		if (MPI_SUCCESS != MPI_Bcast(buffer, INT_MAX, datatype, root, comm))
		{
			DEBUG("Failed to broadcast %d count\n", INT_MAX);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		count -= INT_MAX;
		buffer += INT_MAX;
	}

	if (MPI_SUCCESS != MPI_Bcast(buffer, count, datatype, root, comm))
	{
		DEBUG("Failed to broadcast %lu count\n", count);
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	return MPI_SUCCESS;
}

void mpi_open_psam_file()
{
	char filename[32];

	snprintf(&filename[0], 32, "psam/iter%dproc%d.pSAM", g_iter, g_rank);

	if (NULL != g_fd)
	{
		DEBUG("Global file descriptor was not NULL; cannot continue.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	if (NULL == (g_fd = fopen(filename, "w")))
	{
		DEBUG("Failed to open file: %s\n", filename);
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
}

void mpi_close_psam_file()
{
	if (0 != fclose(g_fd))
	{
		DEBUG("Failed to close file descriptor\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	g_fd = NULL;
}
