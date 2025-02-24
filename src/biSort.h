/**
 * \file biSort.h
 * \brief Regroup bisort function and filter family gene
 * \author Thomas MAGALHAES, Paul PANGANIBAN, Nizar HDADECH
 * \date August 26th 2015
 */

#ifndef BISORT

#define BISORT
#include "parser.h"

typedef struct Cluster
{
int  process_rank;
size_t pnbreads;
size_t family_gene;
size_t family_mgene;
char * gene_name;
char * mgene_name;

struct Read* read;

struct Cluster* next;
} Cluster;


/**\brief	Group all reads in one linked list, sorted by their gene id then by their mate gene id.
 * 			Keep only discordant Reads(condition 1).
 * 			Regroup Reads by cluster(fusion).
 * \param reads Linked list of all Reads that process has
 * \param[out] pCluster Local clusters of process
 * \param preadNumberByChr Reads number by chromosomes
 * \param rank Process rank
 * \param[out] nbCluster Number of Cluster
 * \return total reads number
*/
size_t biSort(Read ***reads,Cluster** pCluster,size_t *preadNumberByChr,int rank,size_t * nbCluster);

/**
 * \brief Find family id for each Reads gene
 * \param buffer_family Buffer with the gene's families file
 * \param[out] pCluster Linked list of clusters
 * \return void
 */
void getGeneFamily(char * buffer_family, Cluster * pCluster);

/**
 * \brief Find family id for each Reads mate gene and delete clusters from a same family fusion
 * \param buffer_family Buffer with the gene's families file
 * \param[out] pCluster Linked list of clusters
 * \return Number of clusters delete
 */
size_t getMGeneFamily(char * buffer_family, Cluster ** pCluster);

#endif
