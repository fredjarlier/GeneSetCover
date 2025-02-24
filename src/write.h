/**
 * \file write.h
 * \date 27 aout 2015
 * \brief Regroups all methods to write the cluster output
 * \author	Thomas Magalhaes
 * 			Paul Panganiban
 */

#ifndef WRITE_H
	#define WRITE_H

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "parser.h"
#include "clusterCommunication.h"
#include "biSort.h"

/**
 * \brief Writes a cluster list given in parameter to a specific file
 *
 * \param rank Rank of the process
 * \param num_proc Total number of jobs on MPI_COMM_WORLD
 * \param clusters Linked-list of clusters
 * \param count_clusters Number of clusters in the linked-list
 * \param chrNames Reference array to the chromosomes names.
 * \param fileName Name of input file.
 * \param header Header of input file. Has been read at the beginning of the program.
 * \param outputDir Where to place the output
 * \param exon_valid_cluster Array of matching percentage for the exon condition per cluster
 *
 * \return void
 */
void writeCluster(int rank, int num_proc, Cluster* clusters, size_t count_clusters, char** chrNames, char* fileName, char* header, char* outputDir, char* exon_valid_cluster);

/**
 * \brief Reads a file part and stocks it in data
 *
 * \param rank Rank of the process
 * \param num_proc Total number of jobs on MPI_COMM_WORLD
 * \param clusters Linked-list of clusters
 * \param[out] data Result string. Will contain the read reads.
 * \param[out] clSize Reference array to the cluster size. Will be set.
 * \param count_clusters Number of clusters in the linked-list
 * \param fileName Name of input file.
 *
 * \return The total size of data
 */
size_t writeCluster_read(int rank, int num_proc, Cluster* clusters, char **data, size_t *clSize, size_t count_clusters, char* fileName);

/**
 * \brief Apply a max to the reads forward coord and a min to the reads backward coord
 * This is a simple loop through algorithm.
 *
 * \param[out] forward Array in which the forwards maxs will be keeped
 * \param[out] backward Array in which the backwards mins will be keeped
 * \param c Linked-list of clusters
 * \param count_cluster Number of clusters in the linked-list
 *
 * \return void
 */
void getFusionsBoundaries(size_t* forward, size_t* backward, Cluster* c, size_t count_cluster);

#endif
