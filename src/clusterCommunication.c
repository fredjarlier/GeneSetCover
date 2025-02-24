/**
 * \file clusterCommunication.c
 * \author  Paul PANGANIBAN,Thomas MAGALHAES, Nizar HDADECH
 * \date August 26th 2015
 */

#include "clusterCommunication.h"
#include "mergeSort.h"
#define SIZE_OF_BCAST 100

void parser_ClusterSplit(Cluster * cluster,ClusterSplit** pClusterSplits,size_t pnbCluster,int rank)
{
	ClusterSplit *clustersplit=NULL;
	ClusterSplit* first=NULL;

	if(pnbCluster != 0){
		clustersplit=(ClusterSplit*)malloc(sizeof(ClusterSplit));
		first=clustersplit;
		int i=0;
		for(i=0;i<pnbCluster;i++)
		{
			clustersplit->gene=cluster->read->gene;
			clustersplit->mgene=cluster->read->mgene;
			clustersplit->nbread=cluster->pnbreads;
			clustersplit->maxrank=rank;
			clustersplit->maxnbread=cluster->pnbreads;

			if(i < (pnbCluster-1)){
				clustersplit->next =(ClusterSplit*)malloc(sizeof(ClusterSplit));
			}
			else{
				clustersplit->next=NULL;
			}
			clustersplit=clustersplit->next;
			cluster=cluster->next;
		}
	}

	clustersplit=first;
	*pClusterSplits=clustersplit;
}

void fork_cluster(ClusterSplit * clusterSplit_to_cut, ClusterSplit ** new_clusterSplit_to_add){
	size_t nb_read_max = DEFAULT_MAX_READ_SIZE/DEFAULT_READ_SIZE;

	ClusterSplit * curr = *new_clusterSplit_to_add;

	while(clusterSplit_to_cut->nbread > nb_read_max){
		if(curr == NULL){
			*new_clusterSplit_to_add =(ClusterSplit*)malloc(sizeof(ClusterSplit));
			curr = *new_clusterSplit_to_add;
		}
		else{
			curr -> next = (ClusterSplit*)malloc(sizeof(ClusterSplit));
			curr = curr->next;
		}

		//we copy attributes
		curr->gene = clusterSplit_to_cut->gene;
		curr->mgene = clusterSplit_to_cut->mgene;
		curr->maxnbread = clusterSplit_to_cut->maxnbread;

		//we give him the left over
		if((clusterSplit_to_cut->nbread % (nb_read_max)) > 0){
			curr->nbread= (clusterSplit_to_cut->nbread % (nb_read_max));
			clusterSplit_to_cut->nbread -=curr->nbread;
		}
		else{
			curr->nbread= nb_read_max;
			clusterSplit_to_cut->nbread -=nb_read_max;
		}

		curr->maxrank = clusterSplit_to_cut->maxrank;
		curr->next = NULL;
		printf("new_clusterSplit_to_add->nbread %zu\n",curr->nbread);
		printf("***********************\n");
	}
}


void split_Cluster(ClusterSplit ** pClusterSplits,int numproc)
{

	ClusterSplit* cs_curr=*pClusterSplits;
	ClusterSplit* cs_first = cs_curr;
	size_t nbreads_by_proc[numproc][2];//know the number of read and cluster for each proc
	size_t nbreads_average = 0;
	size_t range = 5; //percent of difference between nb_read by proce and the average
	size_t nb_read_max = DEFAULT_MAX_READ_SIZE/DEFAULT_READ_SIZE;

	int i=0;
	for(i=0;i<numproc;i++){
		nbreads_by_proc[i][0]=0; //nbreads
		nbreads_by_proc[i][1]=0; //nbcluster
	}

	//we need to sort by nb read our clustersplit list by nbreads
	ClusterSplit * cs_fake = (ClusterSplit*)malloc(sizeof(ClusterSplit));
	cs_fake->gene=0;
	cs_fake->mgene=0;
	cs_fake->maxrank=0;
	cs_fake->nbread=0;
	cs_fake->maxnbread=0;
	cs_fake->next=*pClusterSplits;

	cs_curr = cs_fake;
	size_t nb_cs = 0;
	while(cs_curr){
		nb_cs++;
		cs_curr=cs_curr->next;
	}
	cs_curr=cs_fake;

	//we will sort by nbreads (+ to -)
	mergeSortClusterSplit_inverted(cs_curr,nb_cs-1);
	cs_curr= cs_curr->next;
	cs_first = cs_curr;


	//IF THERE ARE CLUSTERS BIGGER THAN THE PROCESS CAPACITY
	//1-we will create smaller clusters from cluster which are bigger than max_nb_read
	ClusterSplit * tmp = cs_curr;
	ClusterSplit* cs_fork=NULL;//use to save the smaller clusters from bigger
	ClusterSplit* cs_fork_first;
	if(cs_curr->nbread > nb_read_max){
		printf("BIG CLUSTER !!!\n");
		int bool_big_cluster = 0;

		while(!bool_big_cluster){
			fork_cluster(tmp,&cs_fork); //Use to create smaller clusters
			tmp=tmp->next;
			if(tmp->nbread <= nb_read_max){
				bool_big_cluster=1;
			}
		}
		cs_fork_first=cs_fork;
		while(cs_fork->next){
			cs_fork=cs_fork->next;
		}
	}
	//2-we add the new clustersplit to the current
	if(cs_fork != NULL){
		printf("ADD SMALL CLUSTER\n");
		cs_fork->next=cs_curr;
		cs_curr=cs_fork_first;

		//WE WILL SORT AGAIN BY NB READ (+ to -) BECAUSE OF THE CS FORK
		nb_cs=0;
		while(cs_curr){
			nb_cs++;
			cs_curr=cs_curr->next;
		}
		cs_fake->next = cs_curr;
		cs_curr=cs_fake;
		mergeSortClusterSplit_inverted(cs_curr,nb_cs);

		cs_curr = cs_curr->next;
		cs_first = cs_curr;
	}//if(cs_fork!=NULL)


	size_t num_proc_local = numproc;
	//We will know how many reads will recv each proc with the current maxrank
	while(cs_curr)
	{
		//if this cluster is full we give this cs to an other process
		if(nbreads_by_proc[cs_curr->maxrank][0] > (nb_read_max-1)){
			cs_curr->maxrank = (cs_curr->maxrank+1)%numproc;
		}
		//else we add it to our nbreads_by_proc
		else{
			nbreads_by_proc[cs_curr->maxrank][0]+=cs_curr->nbread; //nbread is the total number of read
			nbreads_by_proc[cs_curr->maxrank][1]++; //nb cluster

			//the average only count the read of process that are not full
			if( nbreads_by_proc[cs_curr->maxrank][0] <= (nb_read_max-1)
					&&	nbreads_by_proc[cs_curr->maxrank][1] > 1){
				nbreads_average+=cs_curr->nbread;
			}
			cs_curr=cs_curr->next;
		}//else
	}//while cscurr

	cs_curr=cs_first;


	for(i=0;i<numproc;i++){
		printf("Proc %d :: %zu reads -> %zu clusters\n",i,nbreads_by_proc[i][0],nbreads_by_proc[i][1]);
		if(nbreads_by_proc[i][0] == nb_read_max && nbreads_by_proc[i][1] ==1){
			num_proc_local--; //Only process that are not full will participate for sharing cluster
		}
	}

	nbreads_average/=num_proc_local;

	//WE DEFINE THE NUMBER MIN AND MAX OF READ PER PROCESS
	size_t average_min = nbreads_average-(nbreads_average*range/100);

	//We will know the rank with the max and the min nbreads
	int rank_max=-1,rank_min=-1;
	size_t max_nbreads =0, min_nbreads=-1;

	for(i=0;i<numproc;i++){
		//if the process is not full
		if( nbreads_by_proc[i][0] <= (nb_read_max-1) ){
			//if this process have the most reads and have more than one cluster
			// it will send some of those clusters to the rank min
			if(nbreads_by_proc[i][0]>max_nbreads && nbreads_by_proc[i][1]>1 ){
				max_nbreads = nbreads_by_proc[i][0];
				rank_max = i;
			}
			if(nbreads_by_proc[i][0]<min_nbreads){
				min_nbreads = nbreads_by_proc[i][0];
				rank_min=i;
			}
		}
	}//for

	//We will sort by nbreads (- to +)
	cs_fake->next = cs_curr;
	cs_curr = cs_fake;
	mergeSortClusterSplit(cs_curr,nb_cs-1); //cs will be sort by nbreads (- to +)
	cs_curr=cs_curr->next;//we delete the cluster fake
	cs_first=cs_curr;

	//We will share the clustersplit with each process by giving the lowest clusters of the biggest cluster
	int bool=0;//use to know if the lowest process reach the minimum average of reads

	//while the process with the minimum nb_reads has not more reads than the average_minimum
	//OR if there is only one process which has more than one cluster
	while(min_nbreads<average_min && num_proc_local>1){
		bool=0;
		cs_curr=cs_first;

		//we manage clusters between MAXRANK and MINRANK
		while(!bool){
			//we find the lowest cluster of rank max
			while(cs_curr->maxrank != rank_max ){
				cs_curr=cs_curr->next;
			}
			cs_curr->maxrank=rank_min;//we give the lowest cluster of rank_max to rank_min

			//we update our data
			nbreads_by_proc[rank_min][0]+=cs_curr->nbread;
			nbreads_by_proc[rank_max][0]-=cs_curr->nbread;
			nbreads_by_proc[rank_min][1]++;
			nbreads_by_proc[rank_max][1]--;

			//if the process have now enougth reads OR rankmax has just one cluster
			if((nbreads_by_proc[rank_min][0]>=average_min )
					|| nbreads_by_proc[rank_max][1]==1){
				bool=1;
			}
			//else we try to give it an other cluster
			else{
				cs_curr=cs_curr->next;
			}
		}//while bool

		//we calculate the new rank_min and rank_max
		rank_max=-1;
		rank_min=-1;
		max_nbreads=0;
		min_nbreads=-1;

		num_proc_local=numproc;//use to know if there still are process we can exchange clusters
		for(i=0;i<numproc;i++){
			//if process is not full and have more than one cluster
			if(  nbreads_by_proc[i][0] <= (nb_read_max-1) ){

				if(nbreads_by_proc[i][0]>max_nbreads && nbreads_by_proc[i][1]>1){
					max_nbreads = nbreads_by_proc[i][0];
					rank_max = i;
				}
				if(nbreads_by_proc[i][0]<min_nbreads){
					min_nbreads = nbreads_by_proc[i][0];
					rank_min=i;
				}
			}

			else{
				num_proc_local--;
			}
		}//for

	}//while min_nbread < min average
	printf("************AFTER SPLIT CLUSTER*************\n");
	for(i=0;i<numproc;i++){
		printf("\tProc %d :: %zu reads -> %zu clusters\n",i,nbreads_by_proc[i][0],nbreads_by_proc[i][1]);
	}

	cs_fake->next = cs_first;
	cs_curr = cs_fake;
	mergeSortClusterSplit_by_gene(cs_curr,nb_cs-1); //cs will be sort by nbreads (- to +)
	cs_curr=cs_curr->next;//we delete the cluster fake
	cs_first=cs_curr;
	free(cs_fake);
}

size_t delete_clustersplit_not_enougth_reads (ClusterSplit ** pClusterSplits,size_t nb_reads_min){

	ClusterSplit * cs_curr = *pClusterSplits; //use to browse the list of clustersplit
	ClusterSplit * before = NULL;
	ClusterSplit ** first_cluster = pClusterSplits;
	size_t nb_cs_to_supp = 0;//count nb cs we will supp
	size_t nb_cs_to_keep = 0;//count nb cs we will keep

	// RANK 0 will count how many cluster we will keep and delete clusterSplit with nbread =1
	while(cs_curr != NULL){
		if(cs_curr->nbread > nb_reads_min){
			nb_cs_to_keep++;

			before = cs_curr;
			cs_curr = cs_curr->next;
		}

		//we will delete it from our list
		else{
			nb_cs_to_supp++;
			//if it s the first cluster
			if(before == NULL){
				*first_cluster=cs_curr->next;
				free(cs_curr);
				cs_curr=*first_cluster;
			}
			else{
				before->next=cs_curr->next;
				free(cs_curr);
				cs_curr=before->next;
			}

		}

	}//while

	cs_curr=*first_cluster; // we replace the current to the first cluster

	printf("TOTAL CLUSTER TO SUPP (NBREADS < %zu): %zu / %zu\n",nb_reads_min,nb_cs_to_supp,nb_cs_to_supp+nb_cs_to_keep);
	*pClusterSplits =*first_cluster;
	return(nb_cs_to_keep);

}


void update_cluster_with_csupdate(Cluster ** pCluster,ClusterSplit ** cs_update, size_t * nb_reads_local,size_t * nb_cs_local,size_t nb_read_minimum,int rank){
	Cluster *cl_curr = *pCluster;
	Cluster *cl_before = NULL;
	ClusterSplit * cs_curr_update = *cs_update;
	ClusterSplit * cs_curr_before = NULL;
	ClusterSplit * cs_first = NULL;

	size_t gene=0,mgene=0, nb_reads =0, nb_cs = 0;


	//if we have some clusters
	if(cl_curr){
		//We will supp clusters that are not in the cs_update
		while(cl_curr){
			//if this cluster has more than nb_read_minimum reads we know that we will keep it
			gene = cl_curr->read->gene ;
			mgene = cl_curr->read->mgene ;

			if(cl_curr->pnbreads > nb_read_minimum){

				nb_cs++;
				nb_reads+=cl_curr->pnbreads;
				//we move our clusterSplit pointer to the good clustersplit
				//if its not the next, we dont have the next clusterSplit in ours cluster, so we delete it
				while(cs_curr_update->gene != gene || cs_curr_update->mgene != mgene){

					//if this clustersplit is not for us we delete it
					if(cs_curr_update->maxrank != cl_curr->process_rank){

						//if it's the first
						if(!cs_curr_before){
							*cs_update=(*cs_update)->next;
							free(cs_curr_update);
							cs_curr_update = *cs_update;
						}

						else{
							cs_curr_before->next= cs_curr_update->next;
							free(cs_curr_update);
							cs_curr_update = cs_curr_before->next;
						}
					}

					else{
						cs_curr_before = cs_curr_update;
						cs_curr_update =cs_curr_update->next;
					}
				}//while(we find the good clusterSplit)

				cs_curr_before = cs_curr_update;
				cs_curr_update = cs_curr_update->next;

				cl_before = cl_curr;
				cl_curr = cl_curr->next;
			}

			else{
				//We try to find this cluster in our cs_update
				int bool = 0; //0 : cs_update not found yet / 1: find / 2:to delete


				while(bool==0 && cs_curr_update != NULL ){

					//if we keep this cluster we go out the loop
					if(cs_curr_update->gene == gene && cs_curr_update->mgene == mgene ){
						cs_curr_before = cs_curr_update;
						cs_curr_update = cs_curr_update->next;
						bool=1;
					}

					//if its next we look the next and delete the clusterSplit current
					else if(cs_curr_update->gene < gene ||  (cs_curr_update->gene == gene && cs_curr_update->mgene <mgene) ){
						if(cs_curr_update->maxrank != cl_curr->process_rank){

							//if it's the first
							if(!cs_curr_before){
								*cs_update = cs_curr_update->next;
								free(cs_curr_update);
								cs_curr_update = *cs_update;
							}
							else{
								cs_curr_before->next= cs_curr_update->next;
								free(cs_curr_update);
								cs_curr_update = cs_curr_before->next;
							}
						}
						else{
							cs_curr_before = cs_curr_update;
							cs_curr_update =cs_curr_update->next;
						}

					}//else if

					//if we have to delete it
					else if(cs_curr_update->gene > gene ||  (cs_curr_update->gene == gene && cs_curr_update->mgene >mgene) ){
						bool=2;
					}
				}//while (find the cs_update gene mgene)


				//WHAT DO WE DO WITH THIS CLUSTER

				//if we keep this clustersplit (because we found it) we do the next iteration
				if(bool==1){
					nb_cs++;
					nb_reads+=cl_curr->pnbreads;
					cl_before=cl_curr;
					cl_curr=cl_curr->next;
				}

				//if we didnt find it, we delete this cluster !!!
				else if(bool==0 || bool==2 ){

					//if it s the first cluster
					if(!cl_before){
						*pCluster = cl_curr->next;
						free(cl_curr);
						cl_curr= (*pCluster);
					}
					else{
						cl_before->next= cl_curr->next;
						free(cl_curr);
						cl_curr=cl_before->next;
					}
					//FIN DELETE CLUSTER

				}//else if bool 0 or 2
			}//else nbreads <=nb_read_minimum

		}//while we checked all clusters


		//if there are still clustersplit , we delete them

		if(cs_curr_update ){
			if(gene != cs_curr_update->gene || mgene != cs_curr_update->mgene){
				while(cs_curr_update){
					cs_curr_before->next = cs_curr_update->next;
					free(cs_curr_update);
					cs_curr_update = cs_curr_before->next;
				}
			}
		}
	}

	//if we dont have clusters
	//we only keeped clusterSplit which we are maxrank
	else{
		cs_curr_before = NULL;

		while(cs_curr_update){
			//if we delete this clustersplit
			if(cs_curr_update->maxrank != rank){
				if(cs_curr_before){
					cs_curr_before->next = cs_curr_update->next;
					free(cs_curr_update);
					cs_curr_update= cs_curr_before->next;
				}
				else{
					ClusterSplit * cs_tmp = cs_curr_update;
					cs_curr_update=cs_curr_update->next;
					free(cs_tmp);
				}
			}

			//if we keep this cluster
			else{
				if(!cs_curr_before){
					cs_first=cs_curr_update;
					cs_curr_before=cs_curr_update;
					cs_curr_update=cs_curr_update->next;
				}
				else{
					cs_curr_before=cs_curr_update;
					cs_curr_update=cs_curr_update->next;
				}
				nb_cs++;
			}

		}//while cs_curr_update

		(*cs_update) = cs_first;

	}




	//printf("TOTAL NB READS : %zu // TOTAL NB CLUSTERS %zu\n",nb_reads,nb_cs);

	*nb_reads_local=nb_reads;
	*nb_cs_local=nb_cs;

}


void same_clusterSplit_all(ClusterSplit ** cs_update,size_t *nb_cs_keeped ,int rank, int num_proc)
{
	ClusterSplit *cs_curr = *cs_update; //use to browse our list of clusterSplit
	size_t nb_cs= (*nb_cs_keeped);


	//Création DATATYPE dT_CL_SPLIT
	int count = 3;
	int blocklens[count];
	MPI_Aint indices[count];
	MPI_Datatype old_types[count];
	MPI_Datatype dt_clsplit;
	indices[0]=0;
	old_types[0]=MPI_UNSIGNED_LONG;
	blocklens[0]=4;

	MPI_Type_extent(MPI_UNSIGNED_LONG, &indices[1]);
	indices[1]=indices[1]*blocklens[0];
	old_types[1]=MPI_INT;
	blocklens[1]=1;

	indices[2]=sizeof(ClusterSplit);
	old_types[2]=MPI_UB;
	blocklens[2]=1;

	MPI_Type_struct(count, blocklens, indices, old_types, &dt_clsplit);
	MPI_Type_commit(&dt_clsplit);


	//A flag to tell sending is over
	ClusterSplit * end_cs = (ClusterSplit*)malloc(sizeof(ClusterSplit));
	end_cs->gene=0;
	end_cs->mgene=0;
	end_cs->maxnbread=0;
	end_cs->maxrank=0;
	end_cs->nbread=0;



	//we will supp all clusterplit_update for each process which are not 0
	if(rank){
		ClusterSplit * tmp = NULL;
		while(cs_curr && cs_curr->next){
			tmp = cs_curr->next;
			free(cs_curr);
			cs_curr=tmp;
		}
		//we free the last and malloc with the size of the final list of clustersplit
		if(cs_curr)
			free(cs_curr);

		cs_curr=(ClusterSplit*)malloc(sizeof(ClusterSplit));
		cs_curr->next = NULL;
		*cs_update = cs_curr;
	}

	ClusterSplit* bcastBuf = (ClusterSplit*)malloc(sizeof(ClusterSplit)*SIZE_OF_BCAST);
	ClusterSplit* THE_BEFORE=NULL;
	//cs_curr = *cs_update;
	size_t index = 0;

	while(index<nb_cs)
	{
		if(!rank)
		{
			int i =0;

			for(i=0; i<SIZE_OF_BCAST; i++)
			{
				if(cs_curr != NULL)
				{
					bcastBuf[i] = *cs_curr;
					cs_curr = cs_curr->next;
					index++;
				}
				else
				{
					bcastBuf[i]=(*end_cs);
				}
			}
		}

		MPI_Bcast(bcastBuf,SIZE_OF_BCAST,dt_clsplit,0,MPI_COMM_WORLD);

		if(rank){
			int j = 0;
			int bool = 0;

			//we will add the clusterSplit recv to our cluster_update
			while(!bool && j<(SIZE_OF_BCAST)){
				if(bcastBuf[j].gene == 0 && bcastBuf[j].mgene == 0){
					free(cs_curr);
					cs_curr=NULL;
					THE_BEFORE->next = NULL;
					bool=1;
				}

				else{
					*cs_curr=bcastBuf[j];
					if(index<(nb_cs-1)){
						cs_curr->next=(ClusterSplit*)malloc(sizeof(ClusterSplit));
						cs_curr->next->next = NULL;
					}

					else{
						cs_curr->next=NULL;
					}

					THE_BEFORE=cs_curr;
					cs_curr=cs_curr->next;
					index++;
					j++;

				}
			}

		}//if rank

	}//while(index)

	free(bcastBuf);
	free(end_cs);
	MPI_Type_free(&dt_clsplit);
}



size_t affiche_clusterSplit(ClusterSplit * cluster){
	ClusterSplit* c = cluster ;
	printf("************** AFFICHAGE CLUSTERSPLIT ***************\n");
	size_t i=0;
	while(c){
		i++;
		printf(" %zu - %zu /maxrank: %d / maxnbread %zu nbread %zu \n",c->gene,c->mgene,c->maxrank,c->maxnbread,c->nbread);
		c = c ->next;
	}
	printf("number of clusterSplit %zu\n",i);
	return(i);
}


void all_cluster( ClusterSplit ** clusters, int nbProcess , int rank ){
	int levelMax = level(nbProcess-1);
	size2_t childs ;
	ClusterSplit * clusterSplit = *clusters;
	ClusterSplit * clusterSplit_send = NULL;

	ClusterSplit * new =NULL; //ClusterSplit that we recv and add to ours local clusterSplit


	int i;

	//CREATING DATATYPE dT_CL_SPLIT
	int count = 3;
	int blocklens[count];
	MPI_Aint indices[count];
	MPI_Datatype old_types[count];
	MPI_Datatype dt_clsplit;
	indices[0]=0;
	old_types[0]=MPI_UNSIGNED_LONG;
	blocklens[0]=4;

	MPI_Type_extent(MPI_UNSIGNED_LONG, &indices[1]);
	indices[1]=indices[1]*blocklens[0];
	old_types[1]=MPI_INT;
	blocklens[1]=1;

	indices[2]=sizeof(ClusterSplit);
	old_types[2]=MPI_UB;
	blocklens[2]=1;

	MPI_Type_struct(count, blocklens, indices, old_types, &dt_clsplit);
	MPI_Type_commit(&dt_clsplit);

	//For multiple send
	int cs_buffer_size=50, multiple_loop_i = 0; //50 seems to be the best buffer size
	size2_t cs_buffer_recv_size = {0,0}, count_recv = {0,0};
	ClusterSplit *csBuffer=NULL;
	ClusterSplit *csBufferB=NULL;
	ClusterSplit *csBufferS=NULL; //clusterSpplit buffer to send

	ClusterSplit * end_cs = (ClusterSplit*)malloc(sizeof(ClusterSplit)); //A flag to tell sending is over
	end_cs->gene=0;
	end_cs->mgene=0;
	end_cs->maxnbread=0;
	end_cs->maxrank=0;
	end_cs->nbread=0;

	char  a,b; //boolean uses to know if child a or b have finish to send
	a=0; b=0;
	//We cover the tree level by level
	for(i = levelMax  ; i>0 ; i--){

		//if it is a child
		//We will send ours clusterSplit updated by children's clustersplit or not (if it s a leaf)
		if (level(rank) == i){
			clusterSplit_send = clusterSplit;
			csBufferS = (ClusterSplit*) malloc(sizeof(ClusterSplit)*cs_buffer_size);

			//we count buffer size and send it to father process
			while(clusterSplit_send != NULL){
				cs_buffer_recv_size.a++;
				clusterSplit_send = clusterSplit_send->next;
			}
			MPI_Send(&cs_buffer_recv_size.a, 1, MPI_UNSIGNED_LONG, father(rank), 0, MPI_COMM_WORLD);
			clusterSplit_send = clusterSplit;//go back to first clusterSplit

			//we will send ours clusterSplit
			if(clusterSplit_send){
				while(clusterSplit_send != NULL){
					//Multiple send
					for(multiple_loop_i = 0; multiple_loop_i < cs_buffer_size; multiple_loop_i++)
					{
						if(clusterSplit_send!=NULL)
						{
							csBufferS[multiple_loop_i] = *clusterSplit_send;
							clusterSplit_send = clusterSplit_send->next;
						}
						//if we send all we send flag end
						//
						else
						{
							csBufferS[multiple_loop_i] = *end_cs;
							a=1;
							break;
						}
					}
					MPI_Send(csBufferS,cs_buffer_size,dt_clsplit,father(rank),0,MPI_COMM_WORLD);
				}//while
			}
			free(csBufferS);

		}//END if it is a child


		//If it is a father
		//We will recv clusterSplits from ours childrens
		else if( level(rank) == (i-1) ) {

			//define the child of process rank
			child(rank,&childs, nbProcess);

			if(childs.a==0){
				a = 0;
			}
			else {
				a = 1;
			}

			if(childs.b==0){
				b = 0;
			}

			else {
				b = 1;
			}

			//Recv ClusterSplit buffer from ours childs
			if(a){
				MPI_Recv(&cs_buffer_recv_size.a, 1, MPI_UNSIGNED_LONG, childs.a, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if(b){
				MPI_Recv(&cs_buffer_recv_size.b, 1, MPI_UNSIGNED_LONG, childs.b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}


			csBuffer = (ClusterSplit*) malloc(sizeof(ClusterSplit)*( cs_buffer_recv_size.a + (cs_buffer_size-( cs_buffer_recv_size.a % cs_buffer_size)))); //Complete to a cs_buffer_size multiple
			csBufferB = (ClusterSplit*) malloc(sizeof(ClusterSplit)*( cs_buffer_recv_size.b + (cs_buffer_size-( cs_buffer_recv_size.b % cs_buffer_size)))); //Complete to a cs_buffer_size multiple


			//will stop when  child a and b send all their clusters
			while (a!=0 || b!=0){

				if (a!=0 && cs_buffer_recv_size.a!=0){
					//We will recv a clusterSplit from child.a
					MPI_Recv(&(csBuffer[count_recv.a]), cs_buffer_size,dt_clsplit, (int) childs.a , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
					count_recv.a+=cs_buffer_size;
					if(count_recv.a>=cs_buffer_recv_size.a)
					{
						count_recv.a = cs_buffer_recv_size.a;
						a=0;
					}
				}//if a!=0
				else{
					a=0;
				}

				if (b!=0 && cs_buffer_recv_size.b !=0){
					MPI_Recv(&(csBufferB[count_recv.b]), cs_buffer_size,dt_clsplit, (int) childs.b , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
					count_recv.b+=cs_buffer_size;
					if(count_recv.b>=cs_buffer_recv_size.b)
					{
						count_recv.b = cs_buffer_recv_size.b;
						b=0;
					}
				}//if b!=0
				else{
					b=0;
				}

			}//while
			//END RECV ALL INFORMATION ABOUT CLUSTER CHILD


			//Compute everything
			ClusterSplit *cs_list=NULL, *before=NULL;

			//We add these buffers (child a and b) in one listed link (cs_list)
			for(multiple_loop_i = 0; multiple_loop_i < cs_buffer_recv_size.b; multiple_loop_i++)
			{
				new = (ClusterSplit*)malloc(sizeof(ClusterSplit));
				if(before)
					before->next = new;
				if(cs_list == NULL)
					cs_list = new;
				before = new;
				new->gene=csBufferB[multiple_loop_i].gene;
				new->mgene=csBufferB[multiple_loop_i].mgene;
				new->maxnbread=csBufferB[multiple_loop_i].maxnbread;
				new->maxrank=csBufferB[multiple_loop_i].maxrank;
				new->nbread=csBufferB[multiple_loop_i].nbread;
				new->next = NULL;
			}

			for(multiple_loop_i = 0; multiple_loop_i < cs_buffer_recv_size.a; multiple_loop_i++)
			{
				new = (ClusterSplit*)malloc(sizeof(ClusterSplit));
				if(before)
					before->next = new;
				if(cs_list == NULL)
					cs_list = new;
				before = new;
				new->gene=csBuffer[multiple_loop_i].gene;
				new->mgene=csBuffer[multiple_loop_i].mgene;
				new->maxnbread=csBuffer[multiple_loop_i].maxnbread;
				new->maxrank=csBuffer[multiple_loop_i].maxrank;
				new->nbread=csBuffer[multiple_loop_i].nbread;
				new->next = NULL;
			}

			//If we recvd some clusters we'll add it to the list
			if(cs_list){
				before = (ClusterSplit*)malloc(sizeof(ClusterSplit));
				*before = *end_cs;
				before->next = cs_list;
				cs_list = before;
				mergeSortClusterSplit_by_gene(cs_list, cs_buffer_recv_size.a + cs_buffer_recv_size.b); //we do a bi_sort on clusterSplit rcv to make easier updating phase
				before = cs_list;
				cs_list = cs_list->next;
				free(before);

				ClusterSplit* cs_list_curr = cs_list; //clusterSplit current to ADD (clusterSplit recv)
				ClusterSplit* cs_curr = clusterSplit;//clusterSplit current in LOCAL
				ClusterSplit* before_list = NULL;
				before = NULL;

				//We add them all
				while(cs_list_curr!=NULL && cs_curr!=NULL){
					//If we had it yet, we update our clusterSplit
					if(cs_list_curr->gene == cs_curr->gene && cs_list_curr->mgene == cs_curr->mgene)
					{
						cs_curr->nbread+=cs_list_curr->nbread;

						//if new has a rank with more reads than c
						if(cs_list_curr->maxnbread > cs_curr->maxnbread){
							cs_curr->maxnbread = cs_list_curr->maxnbread;
							cs_curr->maxrank = cs_list_curr->maxrank;
						}

						before_list = cs_list_curr;
						cs_list_curr = cs_list_curr->next;
						free(before_list);
					}

					//Else we add it
					else{
						if(cs_list_curr->gene > cs_curr->gene || ( (cs_list_curr->gene==cs_curr->gene) && (cs_list_curr->mgene > cs_curr->mgene) ) ){
							//That cluster is after
							before = cs_curr;
							cs_curr = cs_curr->next;
						}
						else if (cs_list_curr->gene < cs_curr->gene || ( (cs_list_curr->gene==cs_curr->gene) && (cs_list_curr->mgene < cs_curr->mgene) ) ){
							//Backup next
							before_list = cs_list_curr->next;

							//That cluster must be added !
							if(before != NULL){
								before->next = cs_list_curr;
								cs_list_curr->next = cs_curr;
								cs_curr = cs_list_curr;
							}
							else{
								clusterSplit = cs_list_curr;
								clusterSplit->next = cs_curr;
								cs_curr = cs_list_curr;
							}

							//go to the next
							cs_list_curr = before_list;
						}
					}
				} //while

				//If we have new clusters with gene-mgene bigger than local clusters we add them all at the end
				while(cs_list_curr){
					before->next = cs_list_curr;
					before = before->next;
					cs_list_curr = cs_list_curr->next;
				}

		/*		//We go to the end of clusterSplit local
				while(cs_curr){
					before = cs_curr;
					cs_curr = cs_curr->next;
				}


				before = clusterSplit;
				cs_curr = clusterSplit->next;
				while(cs_curr)
				{
					before_list = cs_curr->next;
					if(before->gene == cs_curr->gene && before->mgene == cs_curr->mgene)
					{
						before->nbread+=cs_curr->nbread;
						if(before->maxnbread < cs_curr->maxnbread){
							before->maxnbread = cs_curr->maxnbread;
							before->maxrank = cs_curr->maxrank;
						}
						before->next = cs_curr->next;
						free(cs_curr);
					}
					else
					{
						before = cs_curr;
					}
					cs_curr = before_list;
				}*/
			}

			if(csBuffer != NULL)
			{
				free(csBuffer);
				csBuffer = NULL;
			}
			if(csBufferB != NULL)
			{
				free(csBufferB);
				csBufferB = NULL;
			}
		}//FIN else if it's a father

		cs_buffer_recv_size.a = 0;
		cs_buffer_recv_size.b = 0;
		count_recv.a = 0;
		count_recv.b = 0;

		//We wait the end of send/recv cluster between father and children
		MPI_Barrier( MPI_COMM_WORLD );

	}//for cover the tree

	if(csBuffer != NULL)
	{
		free(csBuffer);
		csBuffer=NULL;
	}
	free(end_cs);
	MPI_Type_free(&dt_clsplit);
	*clusters=clusterSplit;
}


void participate_set_cover(Cluster ** pCluster,Cluster *** cl_participate_gene,Cluster *** cl_participate_mgene, size_t gene_max ,size_t mgene_max,int * bool_participate,size_t nb_minimum_reads, int rank){

	Cluster ** tab_cl_participate_gene = NULL;
	Cluster ** tab_cl_participate_mgene = NULL;
	size_t nb_cl_gene = 0 ;
	size_t nb_cl_mgene = 0;

	Cluster * cl_curr = (*pCluster);

	//we will count how many clusters will participate
	while(cl_curr){

		if(cl_curr->pnbreads>nb_minimum_reads){
			//know if cl_curr has the gene max
			if( (cl_curr->read->gene == gene_max && cl_curr->read->mgene != mgene_max)
					|| (cl_curr->read->mgene == gene_max && cl_curr->read->gene == mgene_max) ){
				nb_cl_gene++;
			}

			//know if cl_curr has the mgene max
			if( (cl_curr->read->gene == mgene_max && cl_curr->read->mgene != gene_max)
					|| (cl_curr->read->mgene == mgene_max && cl_curr->read->gene != gene_max) ){
				nb_cl_mgene++;
			}
		}
		cl_curr=cl_curr->next;
	}
	cl_curr = (*pCluster);

	//if we participate
	if(nb_cl_gene+nb_cl_mgene >0){
		*bool_participate = 1;

		//if we have cluster for gene_max
		//we will put them in tab_cl_participate_gene
		if(nb_cl_gene >0){
			cl_curr= (*pCluster);
			tab_cl_participate_gene = (Cluster**)malloc(sizeof(Cluster)*nb_cl_gene);

			int i = 0 ;
			for(i=0;i<nb_cl_gene;i++){
				tab_cl_participate_gene[i]=NULL;
			}
			i=0;

			while(i<nb_cl_gene && cl_curr){
				if((cl_curr->read->gene == gene_max && cl_curr->read->mgene != mgene_max)
						|| (cl_curr->read->mgene == gene_max && cl_curr->read->gene == mgene_max) ){
					tab_cl_participate_gene[i]=cl_curr;
					i++;
				}
				cl_curr=cl_curr->next;
			}

			(*cl_participate_gene)=tab_cl_participate_gene;
			cl_curr=(*pCluster);
		}

		//if we have cluster for mgene_max
		//we will put them in tab_cl_participate_mgene
		if(nb_cl_mgene >0){
			cl_curr=(*pCluster);
			tab_cl_participate_mgene = (Cluster**)malloc(sizeof(Cluster)*nb_cl_mgene);
			int i = 0 ;
			for(i=0;i<nb_cl_gene;i++){
				tab_cl_participate_mgene[i]=NULL;
			}
			i=0;

			while(i<nb_cl_mgene && cl_curr){
				if((cl_curr->read->gene == mgene_max && cl_curr->read->mgene != gene_max)
						|| (cl_curr->read->mgene == mgene_max && cl_curr->read->gene != gene_max) ){
					tab_cl_participate_mgene[i]=cl_curr;
					i++;
				}
				cl_curr=cl_curr->next;
			}
			(*cl_participate_mgene)=tab_cl_participate_mgene;
			cl_curr=(*pCluster);
		}


	}//if we participate

}

void delete_reads_set_cover(Cluster * cl, size_t * reads, size_t nb_minimum_reads){
	Read * r_curr = cl->read;
	Read * r_first = cl->read;
	Read * r_before = NULL;

	int i = 0;
	//STEP 1 - we delete reads by COORD

	//USE A FAKE FOR MERGE SORT by COORD
	Read * fake =(Read*)malloc(sizeof(Read));
	fake->coord = 0;
	fake->mcoord = 0;
	fake->next = cl->read;
	cl->read = fake;
	r_first =cl->read;
	r_curr = r_first;

	mergeSort(cl->read,cl->pnbreads);

	r_curr = r_first;

	//We delete Reads from "reads" in Clusters cl
	while(r_curr && (reads+i) ){

		//if its the same we delete it , we keep only one read to keep the cluster name
		if(r_curr->coord == reads[i]){

			if(cl->pnbreads>1){
				cl->pnbreads--;

				if(r_before){
					r_before->next = r_curr->next;
					free(r_curr);
					r_curr = r_before->next;
				}
				else{
					r_first = r_curr->next;
					free(r_curr);
					r_curr = r_first;
				}
			}
			//else we keep this read for the cluster name
			else{
				r_curr=r_curr->next;
			}

		}
		//we will compare with the next read
		else if (r_curr->coord > reads[i]){
			i++;
		}
		else if (r_curr->coord <reads[i]){
			r_before = r_curr;
			r_curr=r_curr->next;
		}
	}

	//STEP 2 - we delete reads by MCOORD

	if(cl->pnbreads > nb_minimum_reads){
		//we merge by mcoord
		mergeSort_mcoord(cl->read,cl->pnbreads);

		r_curr = r_first;
		r_before=NULL;
		i=0;

		while(r_curr && (reads+i) ){

			//if its the same we delete it
			if(r_curr->mcoord == reads[i]){

				if(cl->pnbreads>1){
					cl->pnbreads--;

					if(r_before){
						r_before->next = r_curr->next;
						free(r_curr);
						r_curr = r_before->next;
					}
					else{
						r_first = r_curr->next;
						free(r_curr);
						r_curr = r_first;
					}
				}
				//else we keep this read for the cluster name
				else{
					r_curr=r_curr->next;
				}

			}

			//we will compare with the next read
			else if (r_curr->mcoord > reads[i]){
				i++;
			}

			else if (r_curr->mcoord <reads[i]){
				r_before = r_curr;
				r_curr=r_curr->next;
			}
		}
	}

	cl->read = cl->read->next;
	free(fake);

}


void setCover(Cluster ** pCluster,int rank,int numproc, size_t nb_minimum_reads){

	Cluster * cl_curr = *pCluster;
	Cluster * cl_before = NULL;
	Cluster * cl_first = *pCluster;

	size_t gene_max , mgene_max; // name of gene and mgene MAX global

	size_t nb_read_max_local = 0;
	size_t nb_read_max_global=6;// set at 6 to enter in the loop while

	Cluster ** cl_participate_gene = NULL;//Clusters which have the genemax in their fusion
	Cluster ** cl_participate_mgene = NULL;//Clusters which have the Mgenemax in their fusion

	Cluster * cl_max = NULL;//cl_max in local
	Cluster * cl_max_before = NULL;//use to put the cl_max in cl_valid cluster if it s a valid cluster

	Cluster * cl_valid_curr = NULL; // where we put valid clusters
	Cluster * cl_valid_before = NULL;
	Cluster * cl_valid_first = NULL;

	size_t * reads_send = NULL; //list of coord and mcoord of reads from the cluster MAX


	size_t send_cluster[2];//send name of cluster max
	size_t max_cluster[2];// buff to recv the name of cluster max

	t_set_cover n_send,n_recv;//data use to know the process with the most read (USE MPI_REDUCE MAXLOC)

	int bool_rank_chosen =  0;//to know if the process is the master rank
	int bool_participate = 0;//to know if the process will participate to the read communication


	//SET COVER
	while(nb_read_max_global >nb_minimum_reads){
		//INIT DATA
		gene_max=0;
		mgene_max=0;

		bool_rank_chosen = 0;
		bool_participate = 0;

		nb_read_max_local=0;
		nb_read_max_global=0;

		cl_max = NULL;
		cl_max_before=NULL;

		cl_curr = cl_first;
		cl_before = NULL;

		cl_participate_gene = NULL;
		cl_participate_mgene = NULL;

		//FIND THE MAX CLUSTER FOR EACH PROCE
		while(cl_curr){
			if(cl_curr->pnbreads > nb_read_max_local){
				cl_max_before=cl_before;
				cl_max = cl_curr;
				nb_read_max_local = cl_curr->pnbreads;
				gene_max=cl_curr->read->gene;
				mgene_max=cl_curr->read->mgene;
			}
			cl_before = cl_curr;
			cl_curr=cl_curr->next;
		}

		//Prepare MAXLOC
		n_send.nbread_localmax = nb_read_max_local;
		n_send.rank = rank;

		//NOW WE DO A MPI REDUCE MAX LOC TO KNOW WHICH PROCE HAS THE BIGGEST CLUSTER
		MPI_Allreduce(&n_send,&n_recv,1,MPI_LONG_INT,MPI_MAXLOC,MPI_COMM_WORLD);

		nb_read_max_global = n_recv.nbread_localmax;


		if(nb_read_max_global > nb_minimum_reads){
			//IF IT S THE MASTER PROCE
			if(n_recv.rank==rank && nb_read_max_global<DEFAULT_MAX_READ_SIZE){

				bool_rank_chosen = 1;
				bool_participate = 1;

				send_cluster[0]=gene_max;//tell to other process the cluster name
				send_cluster[1]=mgene_max;
			}//if master proce

			//IF ITS NOT THE MASTER PROCE
			else{
				send_cluster[0]=0;
				send_cluster[1]=0;
			}


			//MASTER PROCE TELL WHAT IS THE BIGGEST CLUSTER
			MPI_Allreduce(send_cluster,max_cluster,2,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

			//we want to know if we participate or not
			participate_set_cover(&cl_first,&cl_participate_gene,&cl_participate_mgene,max_cluster[0],max_cluster[1],&bool_participate,nb_minimum_reads,rank);


			int size_comm_split = 0;
			//We will know the size of the cluster communication
			MPI_Allreduce(&bool_participate,&size_comm_split,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

			//if we have to send reads to other process
			if(size_comm_split >0){
				MPI_Comm Comm_participate,Comm_notParticipate;
				int subrank;
				int size_comm = 0;

				//if we participate
				if(bool_participate==1){

					//we set by default the sub rank of the proc chosen to 0
					if(bool_rank_chosen == 1){
						MPI_Comm_split(MPI_COMM_WORLD,bool_participate,0, &Comm_participate);
					}
					else{
						MPI_Comm_split(MPI_COMM_WORLD,bool_participate,1, &Comm_participate);
					}

					MPI_Comm_rank(Comm_participate,&subrank);
					MPI_Comm_size(Comm_participate,&size_comm);

					reads_send =(size_t*)malloc(sizeof(size_t)*nb_read_max_global*2);

					//if this process is the chosen one , it loads their reads
					if(bool_rank_chosen == 1){
						int j = 0;
						Read * r = cl_max->read;
						for(j=0;j<cl_max->pnbreads;j++){
							reads_send[0+j*2]=r->coord;
							reads_send[1+2*j]=r->mcoord;
							r=r->next;
						}
					}

					//we send reads to subprocess
					MPI_Bcast(reads_send,nb_read_max_global*2,MPI_UNSIGNED_LONG,0,Comm_participate);

					//we sort reads recv by coord to make easily the update
					merge_sort_coord(reads_send,nb_read_max_global*2);

					//****DELETE READS IN SUB RANK****
					//we will supp reads from clusters which have genemax or mgenemax
					if(cl_participate_gene){
						int j = 0;
						while(j<(sizeof(cl_participate_gene)/sizeof(size_t))){
							delete_reads_set_cover(cl_participate_gene[j],reads_send,nb_minimum_reads);
							j++;
						}
						free(cl_participate_gene);
					}//if gene

					if(cl_participate_mgene){
						int j = 0;
						while(j<(sizeof(cl_participate_mgene)/sizeof(size_t))){
							delete_reads_set_cover(cl_participate_mgene[j],reads_send,nb_minimum_reads);
							j++;
						}
						free(cl_participate_mgene);
					}//if mgene

					//****END DELETE READS IN SUB RANK****

					//we have done with reads rcv
					free(reads_send);


					//***AD CLUSTER MAX IN VALID CLUSTER***
					//if we are the master process we put the max cluster in valid cluster
					if(bool_rank_chosen == 1){

						cl_valid_curr = cl_max;

						//first we take out cl_max from pcluster
						if(cl_max_before){
							cl_max_before->next = cl_max->next;
						}
						else{
							cl_first = cl_max->next;
						}

						//we add the cl_tmp to the cl_valid cluster
						if(cl_valid_first == NULL){
							//if it's the first valid cluster
							cl_valid_first = cl_valid_curr;
							cl_valid_before = cl_valid_first; //we init the before
							cl_valid_curr->next = NULL;
						}
						//else we have a first and a before
						else{
							cl_valid_before -> next = cl_valid_curr;
							cl_valid_before = cl_valid_curr;
							cl_valid_curr->next = NULL;
						}
					}
					//*** END AD CLUSTER MAX IN VALID CLUSTER***
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Comm_free(&Comm_participate);
				}

				else{
					MPI_Comm_split(MPI_COMM_WORLD,bool_participate, 1, &Comm_notParticipate);
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Comm_free(&Comm_notParticipate);
				}

			}//if size_comm_split >0
		}//if max nbread > nb_minimum_reads
	}


	MPI_Barrier(MPI_COMM_WORLD);

	//if we still have cluster we will delete them
	//we delete the wrong clusters
	size_t count  = 0;
	if(cl_first){
		count = 1;
		Cluster * cl_next = cl_first->next;
		while(cl_next){
			count++;
			Read * r = cl_next->read;
			Read * r_next = r->next;

			while(r_next){
				r->next=r_next->next;
				free(r->next);
				r_next=r->next;
			}
			free(r);

			cl_first->next=cl_next->next;
			free(cl_next);
			cl_next = cl_first->next;
		}
		free(cl_first);
	}

	(*pCluster)=cl_valid_first;
}

