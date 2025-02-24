/**
 * \file exon.h
 * \brief Contains all function to find exon frame
 * \author Paul PANGANIBAN
 * \date August 26th 2015
 */

#ifndef EXON

#define EXON

#include "biSort.h"
#include <sys/stat.h>
#include <stdio.h>
#include "mpi_globals.h"
#include "genes.h"
#include "mergeSort.h"

typedef struct Cluster_chain
{
struct Cluster * cl;
struct Cluster_chain * next;
} Cluster_chain;

/**
 * \brief For each reads of each clusters, it finds forward's and backward's exon frame.
 * 		  Find for each cluster, gene name and mate gene name with their ID
 * \param[out] pCluster Linked list of Clusters
 * \param nb_cluster Number of cluster in "pCluster"
 * \param gene_exon_path Path of File where are exons
 * \param chrNames Tab with all chromosomes name
 * \param rank Process rank
 * \return void
 */
void main_find_exon(Cluster ** pCluster,size_t nb_cluster,int nb_chr, char* gene_exon_path, char ** chrNames,int rank);

/**
 * \brief Get forward's exon frame for each Reads of Clusters clust
 * \param[out] clust One cluster
 * \param start_exon Pointer of the begining of start exon
 * \param end_exon Pointer of the begining of end exon
 * \param exon_frame Pointer of the begining of exon frame
 * \param nb_exon_total Number of exon for this Gene
 * \return void
 */
void get_exon_forward_for_each_reads(Cluster*clust,char* start_exon, char* end_exon, char* exon_frame,int nb_exon_total);

/**
 * \brief Get forward's exon frame for each Reads of each Clusters cl_by_chr
 * \param[out] cl_by_chr Linked list of Clusters_Chain
 * \param buffer_file Char* with the entire exon file
 * \return void
 */
void get_exonForward_for_each_cluster(Cluster_chain * cl_by_chr, char * buffer_file);

/**
 * \brief Get backward's exon frame for each Reads of Clusters clust
 * \param[out] clust One cluster
 * \param start_exon Pointer of the begining of start exon
 * \param end_exon Pointer of the begining of end exon
 * \param exon_frame Pointer of the begining of exon frame
 * \param nb_exon_total Number of exon for this Gene
 * \return void
 */
void get_exon_backward_for_each_reads(Cluster*clust,char* start_exon, char* end_exon, char* exon_frame,int nb_exon_total);

/**
 * \brief Get backward's exon frame for each Reads of each Clusters cl_by_chr
 * \param[out] cl_by_chr Linked list of Clusters_Chain
 * \param buffer_file Char* with the entire exon file
 * \return void
 */
void get_exonBackward_for_each_cluster(Cluster_chain * cl_by_chr, char * buffer_file);
#endif
