/**
 * \file condition2.h
 * \brief Regroup all functions about condition 2
 * \author Paul PANGANIBAN
 * \date August 26th 2015
 */

#ifndef CONDITION2

#define CONDITION2
#include "parser.h"
#include "biSort.h"
#include "mergeSort.h"

/**
 * \brief Calculate abs(c->coord - d->coord)
 * \param c One read
 * \param d One read
 * \return abs(c->coord - d->coord)
 */
size_t dx(Read * c , Read *d);

/**
 * \brief Calculate abs(c->mcoord - d->mcoord)
 * \param c One read
 * \param d One read
 * \return abs(c->mcoord - d->mcoord)
 */
size_t dy(Read* c , Read *d);

/**
 * \brief Delete Clusters that doesnt pass condition 2 with data alpha on pClusters (see Defuse)
 * \param[out] pClusters Linked list of Clusters which we apply condition on them
 * \param alpha Data alpha (see Defuse)
 * \param minimum_number_reads Minimum number of reads in one cluster
 * \return void
 */
void condition2_on_clusters(Cluster ** pClusters , size_t alpha, size_t minimum_number_reads);

/**
 * \brief Delete Reads that doesnt pass condition 2 with data alpha on a linked list of reads (see Defuse)
 * \param[out] read_cluster Linked list of Reads which we apply condition on them
 * \param alpha Data alpha (see Defuse)
 * \return Number of reads delete by condition 2
 */
size_t condtion2_on_read(Read ** reads_cluster, size_t alpha);



#endif
