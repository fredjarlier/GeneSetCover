/**
 * \file read_communication.c
 * \author Paul PANGANIBAN
 * \date August 27th 2015
 */
#include "read_communication.h"
#define MAX_READS_TO_SEND 100 //one read send is equal to 3 size_t (12octets) + 1 unsignedchar (1 octect)  = 13 octects
//On time t to time t+1, one communication between 2 process will be no more than MAX_READS_TO_SEND * 13 octects
#define MAX_READSSEND_TO_SEND 1000
#define MAX_CLUSTERSEND_TO_SEND 1000


static int compare_rate (void const *a, void const *b)
{
	Iteration const *pa = a;
	Iteration const *pb = b;
	if(pa->rate < pb->rate){
		return(-1);
	}
	else{
		return (pa->rate > pb->rate);
	}
}


void know_nbread_send_recv(size_t * nbread_send_by_proc,size_t * nbread_rcv_by_proc,size_t *nbcluster_send_by_proc,size_t *nbcluster_rcv_by_proc, ClusterSplit * cs_update, Cluster * pCluster,int * order_iteration,int rank, int num_proc){
	ClusterSplit * cs_curr = cs_update;
	Cluster * cl_curr = pCluster;
	size_t gene,mgene;
	size_t total_read_send = 0;
	int i = 0 ;

	//FIRST : EACH PROCESS WILL CALCULATE HOW MANY READS THEY WILL SEND TO OHTERS
	for(i=0;i<num_proc;i++){
		nbread_send_by_proc[i]=0;
		nbcluster_send_by_proc[i]=0;
	}

	//We will check for each cluster how many reads we have to send to a proc with the clusterSplit list update
	while(cl_curr){
		gene = cl_curr->read->gene;
		mgene = cl_curr->read->mgene;

		//We find the clusterSplit
		while(cs_curr ->gene != gene || cs_curr->mgene != mgene){
			cs_curr = cs_curr->next;
		}

		//If we didnt find this cluster, there is an error in arguments
		if(cs_curr ->gene != gene || cs_curr->mgene != mgene){
			fprintf(stderr,"ERROR RANK %d ::: ClusterSplit local(%zu-%zu) AND Cluster local(%zu-%zu) are different\n",rank,cs_curr ->gene,cs_curr ->mgene,gene,mgene);
			MPI_Abort(MPI_COMM_WORLD,7);
		}
		//When cs_curr and cl_curr are the same
		//We update the number of read that we will send to max rank
		//if this is our cluster we will send nothing to ourselves
		if(cs_curr->maxrank != rank){
			nbread_send_by_proc[cs_curr->maxrank]+=cl_curr->pnbreads;
			nbcluster_send_by_proc[cs_curr->maxrank]++;
			total_read_send +=cl_curr->pnbreads;
		}
		cl_curr= cl_curr->next;
	}

	//NOW EACH PROCESS WILL KNOW HOW MANY READS IT WILL RCV FROM OTHER PROCESS
	//We will use a SCATTER

	//each process will do a scatter
	for(i=0;i<num_proc;i++){
		MPI_Scatter(nbread_send_by_proc,1,MPI_UNSIGNED_LONG,
				nbread_rcv_by_proc+i,1,MPI_UNSIGNED_LONG,i,MPI_COMM_WORLD);
	}

	for(i=0;i<num_proc;i++){
		MPI_Scatter(nbcluster_send_by_proc,1,MPI_UNSIGNED_LONG,
				nbcluster_rcv_by_proc+i,1,MPI_UNSIGNED_LONG,i,MPI_COMM_WORLD);
	}

	//NOW EACH PROCESS WILL DO NBREAD SEND/NB READ for each process
	// we will use this data to know what iteration is the most stable (send/recv = 1)
	float * rate_send_recv_per_iteration =(float*)malloc(sizeof(float)*(num_proc-1));
	float * new_rate_send_recv_per_iteration =(float*)malloc(sizeof(float)*(num_proc-1));
	int destsend ,srcrecv;

	for(i=0;i<(num_proc-1);i++){
		destsend=(rank+i+1)%num_proc;
		srcrecv=(rank-i-1)%num_proc;
		if(srcrecv<0){
			srcrecv+=num_proc;
		}
		if(nbread_send_by_proc[destsend]==0 || nbread_rcv_by_proc[srcrecv]==0){
			rate_send_recv_per_iteration[i] = 0;
		}
		else{
			rate_send_recv_per_iteration[i] = (float)nbread_send_by_proc[destsend]/(float)nbread_rcv_by_proc[srcrecv];
		}
	}


	//NOW ALL PROCESS WILL DO A ALL REDUCE
	MPI_Allreduce(rate_send_recv_per_iteration,new_rate_send_recv_per_iteration,num_proc-1,
			MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	free(rate_send_recv_per_iteration);

	//we know how stable is one iteration for all process
	Iteration it[num_proc-1];

	for(i=0;i<(num_proc-1);i++){
		new_rate_send_recv_per_iteration[i]/=num_proc;
		new_rate_send_recv_per_iteration[i]-=1;
		if(new_rate_send_recv_per_iteration[i]<0){
			new_rate_send_recv_per_iteration[i] = -new_rate_send_recv_per_iteration[i];
		}
		it[i].rate =new_rate_send_recv_per_iteration[i];
		it[i].num_iteration = i+1;
	}

	free(new_rate_send_recv_per_iteration);

	//We sort it to have the more stable iteration first
	qsort (it, sizeof it / sizeof *it, sizeof *it, compare_rate);

	for(i=0;i<num_proc-1;i++){
		order_iteration[i]=it[i].num_iteration;
	}
}

Cluster_send convert_cluster_to_csend(Cluster * cluster){
	Cluster_send cs;
	cs.gene=cluster->read->gene;
	cs.mgene=cluster->read->mgene;
	cs.nb_read_in_buffer=cluster->pnbreads;
	return cs;
}

Cluster create_cluster_with_clustersend(Cluster_send * cluster_recv , int rank){
	Cluster c ;
	Read * r = (Read*)malloc(sizeof(Read));

	c.next=NULL;
	c.pnbreads=0; //use to know if its a new cluster
	c.process_rank=rank;

	r->gene=cluster_recv->gene;
	r->mgene=cluster_recv->mgene;

	c.read=r;

	return c ;
}

Read_send convert_read_to_rsend(Read * read){
	Read_send rs;
	char c ;
	rs.coord = read->coord;
	rs.mcoord = read->mcoord;

	rs.file_offset = read->offset_source_file;
	rs.offset = read->offset;

	c =0;
	c |=	read->flags.chr;
	c |=	(read->flags.is_mate<<5);
	c |=	(read->flags.replace_gene_with_mgene<<6);
	c |=	(read->flags.left<<7);

	rs.flag=c;

	rs.mchr=read->mchr;

	return rs;
}

Read convert_rsend_to_read(Read_send * rs_send,size_t gene,size_t mgene){

	Read r;

	r.gene=gene;
	r.mgene=mgene;
	r.coord =rs_send->coord;
	r.mcoord=rs_send->mcoord;
	r.offset_source_file = rs_send->file_offset;
	r.offset = rs_send->offset;

	//data that we dont care
	r.next=NULL;
	r.link=NULL;

	r.flags.chr =  ( rs_send->flag ) & (0x1F);
	r.flags.is_mate = ( rs_send->flag  &(0x20) ) >>5;
	r.flags.replace_gene_with_mgene =( rs_send->flag  &(0x40) ) >>6;
	r.flags.left = ( rs_send->flag  &(0x80) ) >>7;

	r.mchr = rs_send->mchr;

	return r;
}


void add_clustersend_and_readsend_to_cluster(Cluster_send * cluster_recv, Read_send * reads_recv,size_t nbclusters_total_recv,Cluster ** pCluster_local,int rank){

	Cluster * cluster_local = *pCluster_local;
	Cluster * cl_curr = cluster_local;
	Cluster *cl_before = NULL;
	Read * r_curr = NULL;
	Read * r_before = NULL;
	int i = 0 ;
	int j = 0;
	int index_read = 0;
	size_t gene=0;
	size_t mgene=0;
	Read * r_to_add = NULL;

	//we will update ours clusters with the news
	for(i=0;i<nbclusters_total_recv;i++){
		cl_curr = *pCluster_local;
		cl_before=NULL;

		//1-we have to find the new cluster in cluster local
		gene = cluster_recv[i].gene;
		mgene = cluster_recv[i].mgene;

		//we compare with the first to have the r_before

		while((cl_curr != NULL) &&
			(cl_curr->read->gene<gene || (cl_curr->read->gene==gene && cl_curr->read->mgene<mgene) )){
			cl_before = cl_curr;
			cl_curr=cl_curr->next;
		}

		//if its a new cluster that we have to add (bigger than we have)
		if(!cl_curr){
			//if we dont have cluster because we send all ours clusters
			//we add it
			if(!cl_before){
				Cluster * cl_tmp =(Cluster*)malloc(sizeof(Cluster));
				*cl_tmp=create_cluster_with_clustersend( cluster_recv+i,rank);
				*pCluster_local = cl_tmp;
				cl_curr=cl_tmp;
			}

			else{
				Cluster * cl_tmp =(Cluster*)malloc(sizeof(Cluster));
				*cl_tmp=create_cluster_with_clustersend( cluster_recv+i,rank);
				cl_before->next=cl_tmp;
				cl_curr=cl_tmp;
			}

		}

		//if this cluster is between two others
		else if(cl_curr->read->gene != gene || cl_curr->read->mgene != mgene ){
			//if its a new first cluster
			if(!cl_before){
				Cluster * cl_tmp =(Cluster*)malloc(sizeof(Cluster));
				*cl_tmp=create_cluster_with_clustersend( cluster_recv+i,rank);
				cl_tmp->next=*pCluster_local;
				*pCluster_local=cl_tmp;
				cl_curr=cl_tmp;
			}
			else{
				Cluster * cl_tmp =(Cluster*)malloc(sizeof(Cluster));
				*cl_tmp=create_cluster_with_clustersend( cluster_recv+i,rank);
				cl_before->next=cl_tmp;
				cl_tmp->next=cl_curr;
				cl_curr=cl_tmp;


			}
		}

		//we will add each read in our cluster cl_cur
		for(j=0;j<(cluster_recv[i].nb_read_in_buffer);j++){
			//we place ours pointers at the begining
			r_curr=cl_curr->read;
			r_before=NULL;

			//we create the read to add
			r_to_add = (Read*)malloc(sizeof(Read));
			(*r_to_add)=convert_rsend_to_read(reads_recv+index_read+j,gene,mgene);


			//if we just create this cluster we replace the first read by this one
			if(cl_curr->pnbreads == 0){
				free(cl_curr->read);
				cl_curr->read = r_to_add;
			}

			//if there are reads, we will place this one in order
			else{
				//FIRST WE WILL COMPARE WITH THE FIRST READ

				//if the read_to_add is smaller than the first we add it here
				if(r_curr && (r_to_add->coord<r_curr->coord || (r_to_add->coord==r_curr->coord && r_to_add->mcoord<r_curr->mcoord )) ){
					r_to_add->next=r_curr;
					cl_curr->read = r_to_add;
				}
				//if its the same read we delete it
				else if (r_to_add->coord==r_curr->coord && r_to_add->mcoord==r_curr->mcoord ){
					//printf("SAME READ !!! %zu - %zu \n",r_to_add->coord,r_to_add->mcoord);
					free(r_to_add);
					cl_curr->pnbreads--; // because we add one after
				}

				//if the read to add his after the first , we have to find where to put it
				else{
					r_before = r_curr;
					r_curr = r_curr->next;
					//we have to find where to put this read
					while(r_curr && (r_to_add->coord>r_curr->coord || (r_to_add->coord==r_curr->coord && r_to_add->mcoord > r_curr->mcoord )) ){
						r_before = r_curr;
						r_curr=r_curr->next;
					}

					//if we go out the loop because of r_curr == NULL then r_to add it s the new last read
					if(!r_curr){
						r_before->next = r_to_add;
					}
					//this read is between the before and the current OR the same as the current
					else{
						//if the current read and the new read are not the same we add this read
						if(r_to_add->coord != r_curr->coord || r_to_add->mcoord != r_curr->mcoord){
							r_before->next = r_to_add;
							r_to_add->next = r_curr;
						}//IF different read

						//if they are the same we delete the new
						else{
							free(r_to_add);
							cl_curr->pnbreads--; // because we add one after
						}
					}//else between to read
				}//else (if the read to add his after the first , we have to find where to put it)

			}//else there are reads

			cl_curr->pnbreads++;
		}//for (add all reads from one cluster)
		index_read+=cluster_recv[i].nb_read_in_buffer;
	}//for cluster
}



void load_read_to_send(Cluster ** pCluster, ClusterSplit** pCs_update,Cluster_send * cs_tosend,Read_send * r_tosend,size_t* nbcluster_send_by_proc, size_t *nbread_send_by_proc,
		size_t * tab_index_cluster,size_t * tab_index_read ,int num_proc, int rank){

	size_t nbread_to_send_total = 0;
	size_t nbcluster_to_send_total = 0;

	int i = 0;
	//we init the tab_index and calculate the total number of reads and clusters to send
	for(i=0;i<num_proc;i++){
		tab_index_cluster[i]=0;
		tab_index_read[i]=0;
		nbread_to_send_total +=nbread_send_by_proc[i];
		nbcluster_to_send_total+=nbcluster_send_by_proc[i];
	}


	//we will set the begining of each process in the buffer
	size_t index_curr_read = nbread_send_by_proc[0];
	size_t index_curr_cluster= nbcluster_send_by_proc[0];

	for(i=1;i<num_proc;i++){

		tab_index_cluster[i]+=index_curr_cluster;
		tab_index_read[i]+=index_curr_read;
		index_curr_cluster+=nbcluster_send_by_proc[i];
		index_curr_read+=nbread_send_by_proc[i];
	}

	//We convert our reads in read_send and our clusters in cluster_send(to use easier MP
	ClusterSplit * cs_curr = *pCs_update;
	ClusterSplit * cs_before = NULL;

	Cluster * cl_curr= *pCluster;
	Cluster * cl_before = NULL;


	//we look all clusters that we will send
	while(cs_curr){
		//if we keep this cluster we skip
		if(cs_curr->maxrank == rank){
			cs_before = cs_curr;
			cs_curr=cs_curr->next;
		}

		//if we have to send this cluster
		else{
			size_t gene = cs_curr->gene;
			size_t mgene = cs_curr->mgene;
			int send_rank = cs_curr->maxrank;
			//we will find the cluster
			while(cl_curr->read->gene != gene || cl_curr->read->mgene != mgene ){
				cl_before = cl_curr;
				cl_curr = cl_curr->next;
			}

			//we convert the cluster to cluster_send in the good index
			cs_tosend[tab_index_cluster[send_rank]]=convert_cluster_to_csend(cl_curr);
			//we update the index
			tab_index_cluster[send_rank]++;


			//Now we will convert all reads from this cluster in our read_to send
			Read * read_curr = cl_curr->read;

			while(read_curr && cl_curr->pnbreads>0){
				r_tosend[tab_index_read[send_rank]]= convert_read_to_rsend(read_curr);

				//we update our index and delete this read
				tab_index_read[send_rank]++;

				cl_curr->read = read_curr->next;
				free(read_curr);
				read_curr = cl_curr->read;
				cl_curr->pnbreads--;
			}

			//we have done with this cluster we can delete it and look the next
			if(cl_before){
				cl_before->next=cl_curr->next;
				free(cl_curr);
				cl_curr=cl_before->next;
			}
			else{
				*pCluster=cl_curr->next;
				free(cl_curr);
				cl_curr=*pCluster;
			}

			//we can delete this cluster split and look the next
			if(cs_before){
				cs_before->next = cs_curr->next;
				free(cs_curr);
				cs_curr=cs_before->next;
			}
			else{
				*pCs_update = cs_curr->next;
				free(cs_curr);
				cs_curr=*pCs_update;
			}

		}//send this cluster

	}//while

	//we replace the index
	for(i=0;i<num_proc;i++){
		tab_index_read[i]-=nbread_send_by_proc[i];
		tab_index_cluster[i]-=nbcluster_send_by_proc[i];
	}

}
void all_send_rcv(Cluster ** pCluster , ClusterSplit ** pCs_update, size_t* nbread_send_by_proc,size_t*nbread_recv_by_proc,size_t* nbcluster_send_by_proc,size_t*nbcluster_recv_by_proc,int * order_iteration,size_t minimum_number_reads,int rank,int num_proc){

	int i = 0;

	size_t nbread_total_send = 0;
	size_t nbread_total_recv = 0;
	size_t nbcluster_total_recv =0;
	size_t nbcluster_total_send = 0;

	//we will create big buffer send and buffer rcv
	//1-we will send two buffer (one for clusters and the other for reads)

	for(i=0;i<num_proc;i++){
		nbread_total_send += nbread_send_by_proc[i];
		nbread_total_recv +=nbread_recv_by_proc[i];
		nbcluster_total_send +=nbcluster_send_by_proc[i];
		nbcluster_total_recv +=nbcluster_recv_by_proc[i];
	}


	//2-WE CREATE THE BUFFERS and the DATATYPES
	Read_send* buff_readsend =(Read_send*)malloc(sizeof(Read_send)*nbread_total_send);
	Cluster_send* buff_clustersend =(Cluster_send*)malloc(sizeof(Cluster_send)*nbcluster_total_send);
	Read_send* buff_readrecv =(Read_send*)malloc(sizeof(Read_send)*nbread_total_recv);
	Cluster_send* buff_clusterrecv =(Cluster_send*)malloc(sizeof(Cluster_send)*nbcluster_total_recv);

	//Création DATATYPE dT_CL_SEND
	int count = 1;
	int blocklens[count];
	MPI_Aint indices[count];
	MPI_Datatype old_types[count];
	MPI_Datatype dt_clsend;
	indices[0]=0;
	old_types[0]=MPI_UNSIGNED_LONG;
	blocklens[0]=3;

	MPI_Type_struct(count, blocklens, indices, old_types, &dt_clsend);
	MPI_Type_commit(&dt_clsend);


	//Création DATATYPE dT_R_SEND
	MPI_Datatype dt_r_send;
	count = 2;
	int blocklens2[count];
	MPI_Aint indices2[count];
	MPI_Datatype old_types2[count];
	indices2[0]=0;
	old_types2[0]=MPI_UNSIGNED_LONG;
	blocklens2[0]=4;

	MPI_Type_extent(MPI_UNSIGNED_LONG, &indices2[1]);
	indices2[1]=indices2[1]*blocklens2[0];
	old_types2[1] = MPI_CHAR;
	blocklens2[1]=2;


	MPI_Type_struct(count, blocklens2, indices2, old_types2, &dt_r_send);
	MPI_Type_commit(&dt_r_send);

	//3-WE WILL LOAD THE SEND_BUFFER
	size_t tab_index_read_send[num_proc];
	size_t tab_index_cluster_send[num_proc];

	load_read_to_send(pCluster,pCs_update,buff_clustersend,buff_readsend,nbcluster_send_by_proc,nbread_send_by_proc,
			tab_index_cluster_send,tab_index_read_send,num_proc,rank);

	//4-WE WILL START THE COMMUNICATION
	int rank_recv_src , rank_send_dst;
	size_t index_recv_read =0; //use to know where to put the read received inside the buffer
	size_t index_recv_cluster=0;
	MPI_Status status;
	//int j =0;


	for(i=0;i<(num_proc-1);i++){
		rank_send_dst = (rank + order_iteration[i])%num_proc;
		rank_recv_src = (rank-order_iteration[i])%num_proc;

		if(rank_recv_src<0){
			rank_recv_src = rank_recv_src+num_proc;
		}


		//we split sending with a limit size of send
		//1-Sending clusters
		int size_send = 0;
		int count_send = (int)nbcluster_send_by_proc[rank_send_dst];
		int index_send = 0;

		int size_recv = 0;
		int count_recv = (int)nbcluster_recv_by_proc[rank_recv_src];
		int index_recv = 0;


		while(count_send>0 || count_recv>0){

			//if we have to send we calculate the size of send
			if(count_send>0){
				if( ((float)count_send/(float)MAX_CLUSTERSEND_TO_SEND)>=1 ){
					size_send = MAX_CLUSTERSEND_TO_SEND;
				}
				else if(((float)count_send/(float)MAX_CLUSTERSEND_TO_SEND)>0){
					size_send = (count_send%MAX_CLUSTERSEND_TO_SEND);
				}
			}


			//if we have to recv we calculate the size of recv
			if(count_recv>0){
				if( ((float)count_recv/(float)MAX_CLUSTERSEND_TO_SEND)>=1 ){
					size_recv = MAX_CLUSTERSEND_TO_SEND;
				}
				else if(((float)count_recv/(float)MAX_CLUSTERSEND_TO_SEND)>0){
					size_recv = (count_recv%MAX_CLUSTERSEND_TO_SEND);
				}
			}


			if(count_send>0 && count_recv>0){
				MPI_Sendrecv(buff_clustersend+tab_index_cluster_send[rank_send_dst]+index_send,size_send,dt_clsend,rank_send_dst,1,
						buff_clusterrecv+index_recv_cluster+index_recv,size_recv,dt_clsend,rank_recv_src,1,MPI_COMM_WORLD,&status);
			}

			//if we just have to send
			else if(count_send>0 && !(count_recv>0)){
				MPI_Send(buff_clustersend+tab_index_cluster_send[rank_send_dst]+index_send,size_send,dt_clsend,rank_send_dst,1,MPI_COMM_WORLD);
			}

			//if we just have to recv
			else if(!(count_send>0) && count_recv>0){
				MPI_Recv(buff_clusterrecv+index_recv_cluster+index_recv,size_recv,dt_clsend,rank_recv_src,1,MPI_COMM_WORLD,&status);
			}
			//we update ours data
			if(count_send>0){
				count_send -=size_send;
				index_send +=size_send;
			}
			if(count_recv>0){
				count_recv -=size_recv;
				index_recv +=size_recv;
			}

		}


		//1-Sending reads
		size_send = 0;
		count_send = (int)nbread_send_by_proc[rank_send_dst];
		index_send = 0;

		size_recv = 0;
		count_recv = (int)nbread_recv_by_proc[rank_recv_src];
		index_recv = 0;


		while(count_send>0 || count_recv>0){

			//if we have to send we calculate the size of send
			if(count_send>0){
				if( ((float)count_send/(float)MAX_READSSEND_TO_SEND)>=1 ){
					size_send = MAX_READSSEND_TO_SEND;
				}
				else if(((float)count_send/(float)MAX_READSSEND_TO_SEND)>0){
					size_send = (count_send%MAX_READSSEND_TO_SEND);
				}
			}


			//if we have to recv we calculate the size of recv
			if(count_recv>0){
				if( ((float)count_recv/(float)MAX_READSSEND_TO_SEND)>=1 ){
					size_recv = MAX_READSSEND_TO_SEND;
				}
				else if(((float)count_recv/(float)MAX_READSSEND_TO_SEND)>0){
					size_recv = (count_recv%MAX_READSSEND_TO_SEND);
				}
			}


			if(count_send>0 && count_recv>0){
				MPI_Sendrecv(buff_readsend+tab_index_read_send[rank_send_dst]+index_send,size_send,dt_r_send,rank_send_dst,1,
						buff_readrecv+index_recv_read+index_recv,size_recv,dt_r_send,rank_recv_src,1,MPI_COMM_WORLD,&status);
			}

			//if we just have to send
			else if(count_send>0 && !(count_recv>0)){
				MPI_Send(buff_readsend+tab_index_read_send[rank_send_dst]+index_send,size_send,dt_r_send,rank_send_dst,1,MPI_COMM_WORLD);
			}

			//if we just have to recv
			else if(!(count_send>0) && count_recv>0){
				MPI_Recv(buff_readrecv+index_recv_read+index_recv,size_recv,dt_r_send,rank_recv_src,1,MPI_COMM_WORLD,&status);
			}
			//we update ours data
			if(count_send>0){
				count_send -=size_send;
				index_send +=size_send;
			}

			if(count_recv>0){
				count_recv -=size_recv;
				index_recv +=size_recv;
			}

		}//wihle reads

		index_recv_cluster+=nbcluster_recv_by_proc[rank_recv_src];
		index_recv_read+=nbread_recv_by_proc[rank_recv_src];
		MPI_Barrier(MPI_COMM_WORLD);
	}


	//we send all read and cluster we can free it
	free(buff_clustersend);
	free(buff_readsend);

	Cluster * pClusters = *pCluster;
	size_t count_read = 0;
	Read * r;

	//Now we add the reads rcv in ours clusters
	add_clustersend_and_readsend_to_cluster(buff_clusterrecv,buff_readrecv,nbcluster_total_recv,pCluster, rank);

	//we can free clusters and reads recv
	free(buff_clusterrecv);
	free(buff_readrecv);


	pClusters = *pCluster;
	Cluster * cl_before = NULL;

	//NOW WE WILL VERIFY THAT ALL CLUSTERS HAVE MORE THAN 5 READS BECAUSE WE DELETE SOME DUPLICATA
	while(pClusters){
		//if this cluster have not enougth read, we delete it
		if(pClusters->pnbreads <= minimum_number_reads){
			//if it's the first
			if(!cl_before){
				//we delete all reads
				if(pClusters->read){
					r = pClusters->read;
					while(r->next){
						Read * tmp = r->next;
						free(r);
						r=tmp;
					}
					free(r);
				}


				Cluster * tmp = pClusters;
				(*pCluster)=pClusters->next;
				free(tmp);
				pClusters=(*pCluster);
			}
			else{
				//we delete all reads
				if(pClusters->read){
					r = pClusters->read;
					while(r->next){
						Read * tmp = r->next;
						free(r);
						r=tmp;
					}
					free(r);
				}

				cl_before->next=pClusters->next;
				free(pClusters);
				pClusters=cl_before->next;
			}
		}

		else{
			r= pClusters->read;
			count_read=0;
			while(r){
				r=r->next;
				count_read++;
			}
			if(count_read != pClusters->pnbreads){
				printf("NOT SAME NUMBER OF READ !!!!!!!!!!!! \n\n");
			}
			cl_before=pClusters;
			pClusters = pClusters->next;

		}
	}//while delete

	MPI_Type_free(&dt_r_send);
	MPI_Type_free(&dt_clsend);
}



