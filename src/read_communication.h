/**
 * \file read_communication.h
 * \brief Contains all functions for grouping reads from a same cluster on the same process
 * \author Paul PANGANIBAN
 * \date August 27th 2015
 */
#ifndef READS_COMM

#define READS_COMM

#include "biSort.h"
#include "clusterCommunication.h"


typedef struct Cluster_information{
	size_t gene;
	size_t mgene;
	size_t nb_read_in_buffer;
	struct Cluster_information * next ;
}Cluster_information;

typedef struct Cluster_send{
	size_t gene;
	size_t mgene;
	size_t nb_read_in_buffer;
}Cluster_send;


typedef struct Read_send{
	size_t coord;
	size_t mcoord;

	size_t offset;
	size_t file_offset;

	char flag;
	char mchr;

}Read_send;


/**
 * \brief Use to prepare sending/receiving "Reads" for each process.
 * 		  We will know the order of the most stable iteration of send/recv between Process
 * 		  We will know for each process how many Reads and Clusters it will send to/rcv from others
 * \param[out] nbread_send_by_proc Array of size_t which will contain how many Reads, process local will send to others
 * \param[out] nbread_rcv_by_proc Array of size_t which will contain how many Reads, process local will receive from others
 * \param[out] nbcluster_send_by_proc Array of size_t which will contain how many Clusters, process local will send to others
 * \param[out] nbcluster_rcv_by_proc Array of size_t which will contain how many Clusters, process local will receive from others
 * \param[out] cs_update Linked list of Clusterplit. Use to know where send a cluster.
 * \param[out] pCluster Linked list of Clusters local (cs_update and pCluster must have same clusters in the same order )
 * \param[out] order_iteration Array organized by most stable Iteration first
 * \param[in] rank Process rank
 * \param[in] num_proc Number of process
 * \return
 */
void know_nbread_send_recv(size_t * nbread_send_by_proc,size_t * nbread_rcv_by_proc,size_t *nbcluster_send_by_proc,size_t *nbcluster_rcv_by_proc, ClusterSplit * cs_update, Cluster * pCluster,int * order_iteration,int rank, int num_proc);

/**
 * \brief Convert a struct Cluster to struct Cluster_send
 * \param cluster Cluster to convert
 * \return cluster converted in Cluster_send
 */
Cluster_send convert_cluster_to_csend(Cluster * cluster);
/**
 * \brief Convert a struct* cluster_recv to struct Cluster
 * \param cluster_recv Cluster_send to convert
 * \param rank Process rank
 * \return cluster_recv converted in struct Cluster
 */
Cluster create_cluster_with_clustersend(Cluster_send * cluster_recv, int rank);

/**
 * \brief Convert a struct* Read to a struct Read
 * \param read to convert
 * \return read converted in struct Read_send
 */
Read_send convert_read_to_rsend(Read * read);

/**
 * \brief Convert astruct struct * Read_send to struct Read
 * \param rs_send to convert
 * \param gene Gene ID for the converted Read_send
 * \param mgene Mate gene ID for the converted Read_send
 * \return rs_send converted in struct Read
 */
Read convert_rsend_to_read(Read_send * rs_send,size_t gene,size_t mgene);


/**
 * \brief Convert a linked list of struct * Cluster_send and their struct* Read_send in strcut Cluster and Read.
 * 		  Add them in linked list of struct* Clusters.
 * \param[in] cluster_recv Linked list of Cluster_send to add in pCluster_local
 * \param[in] reads_recv Linked list of Read_send (reads of cluster_recv) to add in pCluster_local
 * \param[in] nbclusters_total_recv Number of elements in cluster_recv
 * \param[in,out] pCluster_local Linked list of Cluster to update
 * \param[in] rank Process rank
 * \return void
 */
void add_clustersend_and_readsend_to_cluster(Cluster_send * cluster_recv, Read_send * reads_recv,size_t nbclusters_total_recv,Cluster ** pCluster_local,int rank);

/**
 * \brief Load buffers of Cluster_send and Read_send with a list of Cluster and a list of ClusterSplit.
 * 		  Clusters loaded in buffer of Cluster_send are deleted.
 * \param[in,out] pCluster Linked list of Cluster. Clusters which belong to other process than process "rank" will be delete.
 * \param[in,out] pCs_update Linked list of ClusterSplit. ClusterSplits which belong to other process than process "rank" will be delete.
 * \param[out] cs_tosend Linked list of Cluster_send that we will load
 * \param[out] r_tosend Linked list of Read_send that we will load
 * \param[in] nbcluster_send_by_proc Array with the number of Cluster to send to others
 * \param[out] nbread_send_by_proc Array with the number of Read to send to others
 * \param[in,out] tab_index_cluster Array to know where is the first Cluster_send in cs_to_send for each process
 * \param[in,out] tab_index_read Array to know where is the first Read_send in r_to_send for each process
 * \param[in] num_proc Number of process
 * \param[in] rank Process rank
 * \return void
 */
void load_read_to_send(Cluster ** pCluster, ClusterSplit** pCs_update,Cluster_send * cs_tosend,Read_send * r_tosend,size_t* nbcluster_send_by_proc, size_t *nbread_send_by_proc,size_t * tab_index_cluster,size_t * tab_index_read ,int num_proc,int rank);

/**
 * \brief
 * \param[]
 * \param[]
 * \param[]
 * \param[]
 * \param[]
 * \return void
 */
void all_send_rcv(Cluster ** pCluster , ClusterSplit ** pcs_update, size_t* nbread_send_by_proc,size_t*nb_read_recv_by_proc,size_t* nbcluster_send_by_proc,size_t*nbcluster_recv_by_proc,int * order_iteration,size_t minimum_number_reads,int rank,int num_proc);
#endif

