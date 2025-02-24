/**
 * \file clusterCommunication.h
 * \brief Regroup all function for communication of clusters and set_cover
 * \author  Paul PANGANIBAN,Thomas MAGALHAES, Nizar HDADECH
 * \date August 26th 2015
 */
#ifndef CLUSTERCOMMUNICATION
	#define CLUSTERCOMMUNICATION

#include "parser.h"
#include "biSort.h"
#include "tree.h"

#define DEFAULT_READ_SIZE 72 //octets (before382)
#define DEFAULT_MAX_READ_SIZE 1000000000//octets

typedef struct Iteration{
	float rate;
	int num_iteration;
	struct Iteration * next ;
}Iteration;

typedef struct ClusterSplit
{
	size_t gene;
	size_t mgene;
	size_t nbread;
	size_t maxnbread;
	int  maxrank;
	struct ClusterSplit * next; 

}ClusterSplit;

typedef struct t_set_cover
{
	size_t nbread_localmax;
	int rank;
}t_set_cover;

/**
 * \brief Parse Cluster into ClusterSplit structure
 * \param cluster Linked list of cluster that we will parse
 * \param[out] pClusterSplits Parsing "cluster" result
 * \param pnbCluster
 * \param rank
 * \return void
 */
void parser_ClusterSplit(Cluster * cluster,ClusterSplit** pClusterSplits,size_t pnbCluster,int rank);



/**
 * \brief Split a cluster on some process if this cluster is bigger than capacity of one process
 * \param clusterSplit_to_cut ClusterSplit to cut
 * \param[out] new_clusterSplit_to_add Chained list of ClusterSplit from "clusterSplit_to_cut"
 * \return void
 */
void fork_cluster(ClusterSplit * clusterSplit_to_cut, ClusterSplit ** new_clusterSplit_to_add);

/**
 * \brief Define which clusters will belong for each process
 * \param[out] pclusterSplits Linked list with all clusters from all process
 * \param numproc Number of process
 * \return void
 */
void split_Cluster(ClusterSplit ** pClusterSplits,int numproc);

/**
 * \brief Delete ClusterSplit with not enougth reads
 * \param[out] pClusterSplits Linked list of all ClusterSplit of all process
 * \param nb_min_reads Minimum number of reads that one cluster contains
 * \return Number of ClusterSplit deleted
 */
size_t delete_clustersplit_not_enougth_reads (ClusterSplit ** pClusterSplits,size_t nb_min_reads);

/**
 * \brief All process will know all existing clusters (ClusterSplit)
 * \param[out] cs_update Linked list of ClusterSplit (global)
 * \param nb_cs_keeped Number of ClusterSplit keeped
 * \param rank Process rank
 * \param num_proc Number of process
 * \return void
 */
void same_clusterSplit_all(ClusterSplit ** cs_update,size_t *nb_cs_keeped ,int rank,int num_proc);

/**
 * \brief Delete Clusters from pCluster that are not in cs_update
 * \param[out] pCluster Linked list with Cluster that we will compare to "cs_update"
 * \param cs_update Clusters (ClusterSplit) that we will keeped
 * \param[out] nb_reads_local Local number of reads after this function
 * \param[out] nb_cs_local Local number of cluster after this function
 * \param nb_read_minimum Minimum number of read in one cluster
 * \param rank Process rank
 * \return void
 */
void update_cluster_with_csupdate(Cluster ** pCluster,ClusterSplit ** cs_update, size_t * nb_reads_local,size_t * nb_cs_local,size_t nb_read_minimum,int rank);

/**
 * \brief Print all ClusterSplit from a Linked List
 * \param cluster Linked List of ClusterSplit to print
 * \return Number of ClusterSplit
 */
size_t affiche_clusterSplit(ClusterSplit * cluster);

/**
 * \brief Merge all clusterSplit on each process to process 0.
 * \param clusterSplit Linked list of ClusterSplit local
 * \param nbProcess Number of process
 * \param rank Process rank
 * \return void
 */
void all_cluster( ClusterSplit ** clusterSplit, int nbProcess , int rank );

/**
 * \brief Set bool_participate to 1 if we have clusters with gene_max OR mgene_MAX
 * \param pCluster Linked list of Clusters Local
 * \param[out] cl_participate_gene Clusters from pClusters which has gene_max
 * \param[out] cl_participate_mgene Clusters from pClusters which has mgene_max
 * \param gene_max Gene of cluster choosen
 * \param mgene_max Mate gene of cluster choosen
 * \param[out] bool_participate Boolean set to 1 if this process has a cluster with gene_max or mgene_max
 * \param nb_minimum_reads Minimum number of reads
 * \param rank Process rank
 * \return void
 */
void participate_set_cover(Cluster ** pCluster,Cluster *** cl_participate_gene,Cluster *** cl_participate_mgene, size_t gene_max ,size_t mgene_max,int * bool_participate,size_t nb_minimum_reads, int rank);


/**
 * \brief Delete all Reads from reads_to_delete in each Clusters cl
 * \param[out] cl Linked List of Cluster where we will delete reads from reads_to_delete for each one
 * \param reads_to_delete Array of read's coord and read's mcoord that we will delete in each Clusters cl
 * \param nb_minimum_reads Minimum number of reads
 *
 */
void delete_reads_set_cover(Cluster * cl, size_t * reads_to_delete,size_t nb_minimum_reads);


/**
 * \brief Delete Clusters in pCluster, using the maximum parsimony solution (see Defuse documentation)
 * \param[out] pCluster Linked list of clusters local
 * \param rank Process rank
 * \param numproc Number of process
 * \param nb_minimum_reads Minimum number of reads
 * \return void
 */
void setCover(Cluster ** pCluster,int rank,int numproc,size_t nb_minimum_reads);


#endif
