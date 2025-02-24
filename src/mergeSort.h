/**
 * \file mergeSort.h
 * \brief Regroup all functions of sorting
 * \author Nicolas Fedy to line 33. Then Thomas MAGALHAES, Paul PANGANIBAN, Nizar HDADECH
 * \date August 26th 2015
 */

#ifndef MERGESORT_H
	#define MERGESORT_H


#include <stdlib.h>
#include "parser.h"
#include "clusterCommunication.h"


/**
 * \brief Sort linked list  by smallest read->coord first (skip first element)
 * \param[out] c Linked list Read to sort
 * \param n Number of reads in c
 * \return Read*c sorted by smallest read->coord first
 */
Read* mergeSort(Read* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Read
 * \param p Number of elements of c
 * \param[out] d Linked list of Read
 * \param q Number of elements of d
 * \return Read*new sorted by smallest read->coord first
 */
Read* structMerge(Read* c, size_t p, Read* d, size_t q);


/**
 * \brief Sort linked list  by smallest read->mcoord first (skip first element)
 * \param[out] c Linked list Read to sort
 * \param n Number of reads in c
 * \return Read*c sorted by smallest read->mcoord first
 */
Read* mergeSort_mcoord(Read* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Read
 * \param p Number of elements of c
 * \param[out] d Linked list of Read
 * \param q Number of elements of d
 * \return Read*new sorted by smallest read->mcoord first
 */
Read* structMerge_mcoord(Read* c, size_t p, Read* d, size_t q);

/**
 * \brief Sort linked list c by smallest coord first (skip first element)
 * \param[out] c Linked list of Read_Chain
 * \param n Number of Read_Chain in c
 * \return Read_chain * c sorted by smallest coord
 */
Read_chain* mergeSortReadChain(Read_chain* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Read_Chain
 * \param p Number of elements of c
 * \param[out] d Linked list of Read_Chain
 * \param q Number of elements of d
 * \return Read_Chain*new sorted by smallest coord
 */
Read_chain* structMergeReadChain(Read_chain* c, size_t p, Read_chain* d, size_t q);

/**
 * \brief Sort linked list c by smallest genes first (skip first element)
 * \param[out] c Linked list of Read
 * \param n Number of Read in c
 * \return Read*c sorted by smallest genes first
 */
Read* mergeSortGene(Read* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Read
 * \param p Number of elements of c
 * \param[out] d Linked list of Read
 * \param q Number of elements of d
 * \return Read*new sorted by smallest genes first
 */
Read* structMergeGene(Read* c, size_t p, Read* d, size_t q);


/**
 * \brief Sort linked list c by smallest mate genes first (skip first element)
 * \param[out] c Linked list of Read
 * \param n Number of Read in c
 * \return Read*c sorted by smallest mate genes first
 */
Read* mergeSortMgene(Read* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Read
 * \param p Number of elements of c
 * \param[out] d Linked list of Read
 * \param q Number of elements of d
 * \return Read*new sorted by smallest mate genes first
 */
Read* structMergeMgene(Read* c, size_t p, Read* d, size_t q);


/**
 * \brief Sort linked list c by biggest Clusters first (skip first element)
 * \param[out] c Linked list of Cluster
 * \param n Number of Clusters in c
 * \return Cluster * c sorted by biggest Clusters first
 */
Cluster* mergeSortCluster_by_pnbread_inverted(Cluster* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Cluster
 * \param p Number of elements of c
 * \param[out] d Linked list of Cluster
 * \param q Number of elements of d
 * \return Cluster*new sorted by biggest Clusters first
 */
Cluster* structMergeCluster_by_pnbread_inverted(Cluster* c, size_t p, Cluster* d, size_t q);

/**
 * \brief Sort linked list c by smallest mate genes first (skip first element)
 * \param[out] c Linked list of Clusters
 * \param n Number of Clusters in c
 * \return Cluster*c sorted by smallest mate genes first
 */
Cluster* mergeSortCluster_by_Mgene(Cluster* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Cluster
 * \param p Number of elements of c
 * \param[out] d Linked list of Cluster
 * \param q Number of elements of d
 * \return Cluster*new sorted by smallest mate genes first
 */
Cluster* structMergeCluster_by_Mgene(Cluster* c, size_t p, Cluster* d, size_t q);



/**
 * \brief Sort linked list c by smallest genes first (skip first element)
 * \param[out] c Linked list of Clusters
 * \param n Number of Clusters in c
 * \return Cluster*c sorted by smallest genes first
 */
Cluster* mergeSortCluster_by_Gene(Cluster* c, size_t n);

/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Cluster
 * \param p Number of elements of c
 * \param[out] d Linked list of Cluster
 * \param q Number of elements of d
 * \return Cluster*new sorted by smallest genes first
 */
Cluster* structMergeCluster_by_Gene(Cluster* c, size_t p, Cluster* d, size_t q);


/**
 * \brief Sort linked list c by smallest ClusterSplit first (skip first element)
 * \param[out] c Linked list of ClusterSplit
 * \param n Number of ClusterSplit in c
 * \return ClusterSplit*c sorted by smallest ClusterSplit first
 */
ClusterSplit* mergeSortClusterSplit(ClusterSplit* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of ClusterSplit
 * \param p Number of elements of c
 * \param[out] d Linked list of ClusterSplit
 * \param q Number of elements of d
 * \return ClusterSplit*new sorted by smallest ClusterSplit first
 */
ClusterSplit* structMergeClusterSplit(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q);


/**
 * \brief Sort linked list c by biggest ClusterSplit first (skip first element)
 * \param[out] c Linked list of ClusterSplit
 * \param n Number of ClusterSplit in c
 * \return ClusterSplit*c sorted by biggest ClusterSplit first
 */
ClusterSplit* mergeSortClusterSplit_inverted(ClusterSplit* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of ClusterSplit
 * \param p Number of elements of c
 * \param[out] d Linked list of ClusterSplit
 * \param q Number of elements of d
 * \return ClusterSplit*new sorted by biggest ClusterSplit first
 */
ClusterSplit* structMergeClusterSplit_inverted(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q);

/**
 * \brief Sort linked list c by smallest gene first (skip first element)
 * \param[out] c Linked list of ClusterSplit
 * \param n Number of ClusterSplit in c
 * \return ClusterSplit * c sorted by smallest gene first
 */
ClusterSplit* mergeSortClusterSplit_by_gene(ClusterSplit* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of ClusterSplit
 * \param p Number of elements of c
 * \param[out] d Linked list of ClusterSplit
 * \param q Number of elements of d
 * \return ClusterSplit*new sorted by smallest gene first
 */
ClusterSplit* structMergeClusterSplit_by_gene(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q);

/**
 * \brief Sort linked list c by smallest rank first (skip first element)
 * \param[out] c Linked list of ClusterSplit
 * \param n Number of ClusterSplit in c
 * \return ClusterSplit * c sorted by smallest rank first
 */
ClusterSplit* mergeSortClusterSplit_by_rank(ClusterSplit* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of ClusterSplit
 * \param p Number of elements of c
 * \param[out] d Linked list of ClusterSplit
 * \param q Number of elements of d
 * \return ClusterSplit*new sorted by smallest rank first
 */
ClusterSplit* structMergeClusterSplit_by_rank(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q);


/**
 * \brief Sort linked list c by smallest rate first (skip first element)
 * \param[out] c Linked list of Iteration
 * \param n Number of Iteration in c
 * \return Iteration*c sorted by smallest rate first
 */
Iteration* mergeSort_iteration(Iteration* c, size_t n);
/**
 * \brief Merge linked lists c and d.
 * 		  Sort new linked list by
 * \param[out] c Linked list of Iteration
 * \param p Number of elements of c
 * \param[out] d Linked list of Iteration
 * \param q Number of elements of d
 * \return Iteration*new sorted by smallest rate first
 */
Iteration* structMerge_iteration(Iteration* c, size_t p, Iteration* d, size_t q);

/**
 * \brief Use for merge_sort_coord (see MergeSort algorithm)
 */
void merge_coord (size_t *a, size_t n, size_t m);
/**
 * \brief MergeSort a tab of size_t with n elements
 * \param[in,out] a Tab of size_t to sort
 * \param[in] n Number of elements in a
 * \return void
 */
void merge_sort_coord (size_t *a, size_t n);

#endif
