/**
 * \file parallelMergeSort.c
 * \brief Main function
 * \author Frederic JARLIER, Thomas MAGALHAES, Paul PANGANIBAN, Nizar HDADECH.
 * 		   Nicolas FEDY, Leonor SIROTTI for mergeSort part.
 * \date August 27th 2015
 */

/**
 * \mainpage
 *
 * The project name is Big Data processing application, which consists on processing a high amount of information the fastest and more efficient way possible. Our project is biologically oriented, it is about processing genomic type of information. Its goal is to find break points or points of fusion in the cell, which will help scientists to understand the cause of some type of cancer.
The project has been realized within the context of the third degrees "Tutored project" module. The Project has been supervised by Mr. Frederic Jarlier bioinformatics engineer at Curie institute.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <mpi.h>
#include <assert.h>
#include <fcntl.h>
#include <inttypes.h>
#include <libgen.h>
#include <unistd.h>

#include "mpi_globals.h"
#include "write.h"
#include "mergeSort.h"
#include "parser.h"

#include "clusterCommunication.h"
#include "duplicata.h"
#include "genes.h"
#include "read_communication.h"
#include "condition2.h"
#include "exon.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define	PRIoff PRId64

#define	MPI_OFF_T MPI_LONG_LONG_INT

#define DEFAULT_MAX_SIZE 6000000000 //Default capacity for one process: 6G
#define DEFAULT_ALPHA 200
#define DEFAULT_INBUF_SIZE  (512*1024*1024) //4G is the best size for us

int main (int argc, char *argv[]){

	MPI_File mpi_filed;
	FILE* time_output_file;
	int num_proc, rank;
	int nbchr, i;
	int ierr, errorcode = MPI_ERR_OTHER;
	char *file_name, *output_dir, *read_name, *time_dir;
	char *header = NULL;
	char* chr_path, *gene_family_path, *gene_exon_path;
	char **chrNames = NULL;
	struct stat fileSize;
	unsigned int headerSize=0;
	unsigned char threshold = 0;
	char paired = 0;
	size_t unmappedSize = 0;
	size_t *readNumberByChr = NULL;
	Read **reads;
	Read_chain **backward;
	size_t *backwardReadNumber;
	clock_t tic, toc,tac;

	// For clustering - From BiSort
	Cluster *pClusters;
	ClusterSplit *clusterSplit;
	size_t pNumberCluster=0, readNumber=0;


	char *rbuf;
	size_t fsiz, lsiz, loff, *goff;


	MPI_Init(&argc,&argv);
	tac = MPI_Wtime();

	if (argc < 7){
		fprintf(stderr, "Invalid arguments.\nShutting down.\n");
		MPI_Finalize();
		return 0;
	}

	//finds out process rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//finds out number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	if(!rank)fprintf(stderr,"Number of processes : %d\n",num_proc);

	//Init arguments
	file_name = argv[1];
	output_dir = argv[2];
	read_name = argv[3];
//	char *host_name_prefix = argv[4]; Doesn't use but we let it if need
	chr_path = argv[5];
	time_dir = argv[6];
	gene_family_path = argv[7];
	gene_exon_path = argv[8];

	if(rank == 0)fprintf(stderr, "file to read : %s\n", file_name);
	if(rank == 0)fprintf(stderr, "output : %s\n", output_dir);
	if(rank == 0)fprintf(stderr, "readname : %s\n", read_name);
	if(rank == 0)fprintf(stderr, "time : %s\n", time_dir);

	//looking for options
	for(i = 0; i < argc; i++){
		if(argv[i][0] == '-'){
			if(argv[i][1] == 'q'){
				threshold = atoi(argv[i+1]);
				if(!rank)fprintf(stderr, "Reads' quality threshold : %d\n", threshold);
			}

			if(argv[i][1] == 'p'){
				paired = 1;
				if(!rank)fprintf(stderr, "Reads are paired\n");
			}
		}
	}

	//open for append the time file
	time_output_file = fopen(time_dir, "a");
	//time_output_file =stderr;
	//we open the input file
	ierr = MPI_File_open(MPI_COMM_WORLD, file_name,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_filed);
 	if (ierr){
		if (rank == 0) DEBUG("%s: Failed to open file in process 0 %s\n", argv[0], argv[1]);
		MPI_Abort(MPI_COMM_WORLD, errorcode);
		exit(2);
	}

	//then we get the file size
	int input_file_size = stat(file_name, &fileSize);

	if (input_file_size == -1){
		fprintf(stderr,"Failed to find file size\n");
		MPI_Abort(MPI_COMM_WORLD, errorcode);
		exit(0);
	}

	input_file_size = (long long)fileSize.st_size;
	if(!rank)fprintf(stderr, "The size of the file is %zu\n", fileSize.st_size);

	/* Get chunk offset and size */
	fsiz = fileSize.st_size;
	lsiz = fsiz / num_proc;
	loff = rank * lsiz;
	size_t lsiz2 = 150*sizeof(char)*1000; // load only the begining of file to check the read name and header
	if(lsiz2>lsiz){
		lsiz2 = lsiz;
	}

	rbuf = (char*)malloc((lsiz2 + 1)*sizeof(char));

	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***READ AT ALL***\n", rank);
	MPI_File_read_at_all(mpi_filed, loff, rbuf, lsiz2, MPI_CHAR, MPI_STATUS_IGNORE);
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN READ AT ALL***\n", rank,(double)(MPI_Wtime()-tic));


	//FIND HEADERSIZE AND CHRNAMES AND NBCHR
	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***HEADER***\n", rank);
	headerSize=find_header(rbuf,rank,&unmappedSize,&nbchr,&header,&chrNames);
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN HEADER %d***\n", rank,(double)(MPI_Wtime()-tic),headerSize);
	free(rbuf);

	MPI_Status status;

	//We place file offset of each process to the begining of one read's line
	goff=init_goff(mpi_filed,headerSize,fileSize.st_size,num_proc,rank);

	//We calculate the size to read for each process
	lsiz = goff[rank+1]-goff[rank];


	//NOW WE WILL PARSE
	int j=0;
	size_t poffset = goff[rank]; //Current offset in file sam
	reads = (Read**)malloc(nbchr*sizeof(Read));//We allocate a linked list of struct for each Chromosome (last chr = unmapped reads)
	readNumberByChr = (size_t*)malloc(nbchr*sizeof(size_t));//Array with the number of reads found in each chromosome
	Read ** anchor = (Read**)malloc(nbchr*sizeof(Read));//Pointer on the first read of each chromosome

	//Init first read
	for(i = 0; i < nbchr; i++){
		reads[i] = (Read*)malloc(sizeof(Read));
		reads[i]->coord = 0;
		reads[i]->mcoord = 0;
		anchor[i] = reads[i];
		readNumberByChr[i]=0;
	}

	toc = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***FULL PARSER***\n", rank);
	char * local_data = NULL; //Where we load file sam

	//We read the file sam and parse
	while(poffset < goff[rank+1]){
		size_t size_to_read = 0;

		//Reading in multiple times because of MPI_File_read limits
		if( (goff[rank+1]-poffset) < DEFAULT_INBUF_SIZE ){
			size_to_read = goff[rank+1]-poffset;
		}
		else{
			size_to_read = DEFAULT_INBUF_SIZE;
		}

		//we load the buffer
		local_data=(char*)calloc(size_to_read+1,sizeof(char));
		MPI_File_read_at_all(mpi_filed, (MPI_Offset)poffset, local_data, size_to_read, MPI_CHAR, &status);

		//we look where is the last line read for updating next poffset
		size_t offset_last_line = size_to_read-1;
		while(local_data[offset_last_line] != '\n'){
			offset_last_line -- ;
		}
		//If it s the last line of file, we place a last '\n' for the function tokenizer
		if(rank == num_proc-1 && poffset+size_to_read == goff[num_proc]){
			local_data[offset_last_line]='\n';
		}

		//Now we parse Read in local_data
		if(paired){
			parser_paired(local_data, rank, poffset, threshold,nbchr, &readNumberByChr, chrNames, &reads);
		}
		else{
			parser_single(local_data, rank, poffset, threshold,nbchr, &readNumberByChr, chrNames, &reads);
		}

		//we go to the next line
		poffset+=(offset_last_line+1);
		free(local_data);
	}
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN FULL PARSER***\n", rank,(double)(MPI_Wtime()-toc));

	free(goff);
	MPI_File_close(&mpi_filed);
	MPI_Barrier(MPI_COMM_WORLD);

	//We set attribute next of the last read and go back to first read of each chromosome
	for(i = 0; i < nbchr; i++){
		reads[i]->next = NULL;
		reads[i] = anchor[i];
	}
	free(anchor);

	//We count how many reads we found
	size_t nb_reads_total =0,nb_reads_global =0;
	for(j=0;j<nbchr;j++){
		nb_reads_total+=readNumberByChr[j];
	}
	MPI_Allreduce(&nb_reads_total,&nb_reads_global,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
	printf("Number of reads on rank %d = %zu/%zu \n",rank,nb_reads_total,nb_reads_global);


	fprintf(time_output_file, "%d (%.2lf)::::: ***START COUNT ***\n", rank, (double)(MPI_Wtime()-tac));

	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***MERGE SORT***\n", rank);
	//we sort each reads per chromosome by smallest coord first
	for(i = 0; i < nbchr; i++){
		if(reads[i] && reads[i]->next && reads[i]->next->next){
			mergeSort(reads[i], readNumberByChr[i]);
		}
	}
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN MERGE SORT ***\n", rank, (double)(MPI_Wtime()-tic));


	size_t total_readNumberByChr[nbchr];
	MPI_Allreduce(readNumberByChr,total_readNumberByChr,nbchr,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
	//we count the new_nbchr
	for(i=0;i<nbchr && rank==0;i++){
		if(!rank)
			printf("Total reads for %s : %zu\n",chrNames[i],total_readNumberByChr[i]);
	}

	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***GetGenes***\n", rank);

	char * chr_path_temp= NULL;
	MPI_File chrGene_file = NULL;
	struct stat chrGene_file_size;
	size2_t found;

	//Find genes names for each chromosome
	for(i = 0; i < nbchr-1; i++){
		//Change file name to the current chromosome
		if(total_readNumberByChr[i]>0){

			//we set path of file
			chr_path_temp=(char*)calloc(255,sizeof(char));
			strcpy(chr_path_temp,chr_path);
			strcat(chr_path_temp,chrNames[i]);
			strcat(chr_path_temp,".gtf");

			//we open the input file
			ierr = MPI_File_open(MPI_COMM_WORLD, chr_path_temp,  MPI_MODE_RDONLY , MPI_INFO_NULL, &chrGene_file);
			if (ierr){
				if (rank == 0)
					DEBUG("Failed to open file in process 0 %s\n",chr_path_temp);
			}

			else{
				//we take file size
				int chr_file_size = stat(chr_path_temp, &chrGene_file_size);
				if (chr_file_size == -1){
					fprintf(stderr,"Failed to find file size\n");
					MPI_Abort(MPI_COMM_WORLD, errorcode);
					exit(0);
				}
				chr_file_size = (long long)chrGene_file_size.st_size;

				//if we have no gene for this chromosome we delete all reads of this chromosome
				if(chr_file_size == 0){
					Read * next = reads[i]->next;
					int j = 0;
					found.a=0;
					found.b=readNumberByChr[i];
					for(j=0;j<readNumberByChr[i];j++){
						reads[i]->next =next->next;
						free(next);
						next=reads[i]->next;
					}
					reads[i]=reads[i]->next;
					readNumberByChr[i]=0;
				}

				//else we will find gene id for each read
				else{
					char * buffer =(char*)calloc(chr_file_size+1,sizeof(char));
					MPI_Status chr_status;
					MPI_File_read_all(chrGene_file,buffer,chr_file_size,MPI_CHAR,&chr_status);

					found.a=0;//represent the number of reads which maped on gene
					found.b=0;//represent the number of reads which didnt map on gene
					found = getGeneAsNumber(buffer,reads[i]);
					readNumberByChr[i]-=found.b;

					reads[i]=reads[i]->next;//skip next because first read is not a real read

					free(buffer);
				}
				MPI_File_close(&chrGene_file);
			}
			free(chr_path_temp);
		}
	}

	//we count again if there are still reads for each chromosome
	MPI_Allreduce(readNumberByChr,total_readNumberByChr,nbchr,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);


	//Prepare backward
	if(paired){

		backwardReadNumber = (size_t*)calloc(nbchr,sizeof(size_t));
		//Using backward to organize our reads by chromosome of backward
		backward = getBackward(reads, nbchr, backwardReadNumber);

		Read_chain * fake = malloc(sizeof(Read_chain));//Only use on mergeSort because mergeSort skip first element
		//Sort backward by smallest coord first
		for(i = 0; i < nbchr-1; i++){
			if(backward[i] && backward[i]->next && backward[i]->next->next){
				fake->next=backward[i];
				mergeSortReadChain(fake, backwardReadNumber[i]);
				backward[i]=fake->next;
			}
		}
		free(fake);

		//We open each chromosome file to find backward's gene
		for(i=0;i<(nbchr-1);i++){
			if(total_readNumberByChr[i]>0){
				//we set path of file
				chr_path_temp=(char*)calloc(255,sizeof(char));
				strcpy(chr_path_temp,chr_path);
				strcat(chr_path_temp,chrNames[i]);
				strcat(chr_path_temp,".gtf");

				//we open the input file
				ierr = MPI_File_open(MPI_COMM_WORLD, chr_path_temp,  MPI_MODE_RDONLY , MPI_INFO_NULL, &chrGene_file);
				if (ierr){
					if (rank == 0)
						DEBUG("Failed to open file in process 0 %s\n",chr_path_temp);
				}

				else{
					//we take file size
					int chr_file_size = stat(chr_path_temp, &chrGene_file_size);
					if (chr_file_size == -1){
						fprintf(stderr,"Failed to find file size\n");
						MPI_Abort(MPI_COMM_WORLD, errorcode);
						exit(0);
					}
					chr_file_size = (long long)chrGene_file_size.st_size;

					//if we have no gene for this chromosome
					if(chr_file_size == 0){
						found.a=0;
						found.b=readNumberByChr[i];
					}

					//else we will find mate gene id for each backward
					else{
						char * buffer =(char*)calloc(chr_file_size+1,sizeof(char));
						MPI_Status chr_status;
						MPI_File_read_all(chrGene_file,buffer,chr_file_size,MPI_CHAR,&chr_status);

						found = getMGeneAsNumber(buffer,backward[i]);

						free(buffer);
					}

					MPI_File_close(&chrGene_file);
				}
				free(chr_path_temp);
			}
		}

		free(backwardReadNumber);
		assert(backward != NULL);
		freeBackward(backward, nbchr);
	}

	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN GetGenes ***\n", rank, (double)(MPI_Wtime()-tic));

	//#####################################################################################
	//#####################  CLUSTERING START   ###########################################
	//#####################################################################################


	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***BISORT***\n", rank);
	//Regroup our reads in one linked list
	//Group them by fusion
	//Creating linked list of Cluster
	readNumber=biSort(&reads,&pClusters,readNumberByChr,rank,&pNumberCluster);
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN BISORT  nombre de cluster : %zu / nombre de reads %zu***\n", rank, (double)(MPI_Wtime()-tic), pNumberCluster,readNumber);


	//DELETE CLUSTER FROM THE SAME FAMILY
	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***DELETE SAME FAMILY FUSION***\n", rank);

	MPI_File family_file = NULL;
	struct stat family_file_stat;
	size_t nb_same_family = 0 ;

	//We open file with containing Family gene ID
	ierr = MPI_File_open(MPI_COMM_WORLD, gene_family_path,  MPI_MODE_RDONLY , MPI_INFO_NULL, &family_file);
	if (ierr){
		if (rank == 0)
			DEBUG("Failed to open file in process 0 %s\n",gene_family_path);
	}

	else{
		//we take file size
		int family_file_size = stat(gene_family_path, &family_file_stat);
		if (family_file_size == -1){
			fprintf(stderr,"Failed to find file size\n");
			MPI_Abort(MPI_COMM_WORLD, errorcode);
			exit(0);
		}
		family_file_size = (long long)family_file_stat.st_size;

		char * buffer_family_gene =(char*)calloc(family_file_size+1,sizeof(char));
		//We load the entire file in buffer_family_gene
		MPI_File_read_all(family_file,buffer_family_gene,family_file_size,MPI_CHAR,&status);
		//We set gene ID family for each Cluster
		if(pClusters){
			getGeneFamily(buffer_family_gene,pClusters);
		}
		free(buffer_family_gene);

		//WE DO THE SAME FOR MATE GENE ID FAMILY

		//We sort clusters by mate gene to make easier finding id family
		Cluster * fake_cl =(Cluster*)malloc(sizeof(Cluster));//use a fake for mergeSort (skip first element)
		fake_cl->read=(Read*)malloc(sizeof(Read));
		fake_cl->read->gene = 0;
		fake_cl->read->mgene = 0;
		fake_cl->next = pClusters;
		if(fake_cl->next && fake_cl->next->next){
			mergeSortCluster_by_Mgene(fake_cl,pNumberCluster);
		}

		//we skip fake cluster
		pClusters = fake_cl->next;

		//We read again Family gene file because our first before has been edited
		buffer_family_gene =(char*)calloc(family_file_size+1,sizeof(char));
		MPI_File_seek(family_file,0,MPI_SEEK_SET);//go back to the begining of file
		MPI_File_read_all(family_file,buffer_family_gene,family_file_size,MPI_CHAR,&status);

		//we delete clusters with fusion from same family
		if(pClusters){
			nb_same_family=getMGeneFamily(buffer_family_gene,&pClusters);
		}
		pNumberCluster-=nb_same_family;
		free(buffer_family_gene);
		MPI_File_close(&family_file);

		//we sort again ours clusters by gene
		fake_cl->next=pClusters;
		if(fake_cl->next && fake_cl->next->next){
			mergeSortCluster_by_Gene(fake_cl,pNumberCluster);
		}
		pClusters = fake_cl->next;
		free(fake_cl);
	}
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN DELETE SAME FAMILY FUSION ***\n", rank, (double)(MPI_Wtime()-tic));


	tic = MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***DUPLICATA***\n", rank);
	//Delete reads duplicated in local (can have duplicate split on different process)
	size_t nb_duplicata = duplicata(pClusters);
	readNumber-=nb_duplicata;
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN DUPLICATA : %zu reads supp ***\n", rank, (double)(MPI_Wtime()-tic),nb_duplicata);


	//Creation of clusterSplit for all_cluster step
	tic=MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***CREATION CLUSTERSPLIT***\n", rank);
	parser_ClusterSplit(pClusters,&clusterSplit,pNumberCluster,rank);
	fprintf(time_output_file, "%d (%.2lf)::::: *** FIN CREATION CLUSTERSPLIT***\n", rank,(double)(MPI_Wtime()-tic));


	//Tell 0 all existing clusters
	tic=MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***ALLCLUSTER (TREE)***\n", rank);
	all_cluster(&clusterSplit,num_proc,rank);
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN ALLCLUSTER (TREE)***\n", rank,(double)(MPI_Wtime()-tic));


	//Process 0 manages which process will recv which clusters
	tic=MPI_Wtime();
	fprintf(time_output_file, "%d ::::: ***SPLIT CLUSTER (ONLY 0)***\n", rank);
	size_t nb_cs_keeped = 0;
	if(!rank && num_proc>0){
		//Process 0 knows number of read of each cluster, so it deletes all little clussters
		//A fusion must have a minimum number of reads
		nb_cs_keeped=delete_clustersplit_not_enougth_reads(&clusterSplit,5);

		//Process 0 split all clusters equitably between all process
		//We want all reads from a same cluster in one process for the next treatment
		if(nb_cs_keeped>0){
			split_Cluster(&clusterSplit,num_proc);
		}
	}
	fprintf(time_output_file, "%d (%.2lf)::::: ***FIN SPLIT CLUSTER***\n", rank, (double)(MPI_Wtime()-tic) );

	//We will tell every process the size of the new clustersplit list
	MPI_Bcast(&nb_cs_keeped,1,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);

	if(nb_cs_keeped>0){

		//Bcast all Clustersplits (All process will know which clusters to keep)
		tic=MPI_Wtime();
		fprintf(time_output_file, "%d ::::: *** SAME CLUSTERSPLIT FOR ALL***\n", rank);
		Cluster ** ppClusters = &pClusters;
		same_clusterSplit_all(&clusterSplit,&nb_cs_keeped ,rank,num_proc);
		fprintf(time_output_file, "%d (%.2lf)::::: *** FIN SAME CLUSTERSPLIT FOR ALL***\n", rank, (double)(MPI_Wtime()-tic));

		//Update Clusters with the clusterSplit recv from 0
		size_t nb_reads_local=0;
		size_t nb_cs_local=0;
		tic=MPI_Wtime();
		fprintf(time_output_file, "%d ::::: *** UPDATE CLUSTERSPLIT LOCAL WITH NEW ***\n", rank);
		update_cluster_with_csupdate(ppClusters,&clusterSplit,&nb_reads_local,&nb_cs_local,5,rank);
		fprintf(time_output_file, "%d (%.2lf)::::: *** FIN UPDATE CLUSTERSPLIT LOCAL WITH NEW ***\n", rank, (double)(MPI_Wtime()-tic));

		MPI_Barrier(MPI_COMM_WORLD);

		size_t nbcluster_rcv_by_proc[num_proc];
		size_t nbcluster_send_by_proc[num_proc];
		size_t nbread_rcv_by_proc[num_proc];
		size_t nbread_send_by_proc[num_proc];
		int order_iteration[num_proc-1]; //tab use to know what iteration do first for the all_send_recv reads

		//Use to know the number of reads and clusters send between each process
		//And to know what iteration is the most steady to do not have a process overload
		know_nbread_send_recv(nbread_send_by_proc,nbread_rcv_by_proc,nbcluster_send_by_proc,nbcluster_rcv_by_proc,clusterSplit,pClusters,order_iteration,rank,num_proc);

		//We send all cluster to its process max rank
		tic=MPI_Wtime();
		fprintf(time_output_file, "%d ::::: *** ALL SENDRECV READS ***\n", rank);
		all_send_rcv(&pClusters,&clusterSplit,nbread_send_by_proc,nbread_rcv_by_proc,nbcluster_send_by_proc,nbcluster_rcv_by_proc,order_iteration,5,rank,num_proc);
		fprintf(time_output_file, "%d (%.2lf)::::: *** FIN ALL SENDRECV ***\n", rank, (double)(MPI_Wtime()-tic));

		//We apply condition 2 on all clusters (see Defuse)
		tic=MPI_Wtime();
		fprintf(time_output_file, "%d ::::: *** CONDITION 2 ***\n", rank);
		condition2_on_clusters(&pClusters,DEFAULT_ALPHA,5);
		fprintf(time_output_file, "%d (%.2lf)::::: *** FIN CONDITION 2 ***\n", rank, (double)(MPI_Wtime()-tic));

		//We delete fake clusters with set cover algorithm (see Defuse)
		tic=MPI_Wtime();
		fprintf(time_output_file, "%d ::::: *** SET COVER ***\n", rank);
		setCover(&pClusters,rank,num_proc,5);
		fprintf(time_output_file, "%d (%.2lf)::::: *** FIN SET COVER ***\n", rank, (double)(MPI_Wtime()-tic));

		//We print all Valid clusters
		Cluster * test = pClusters;
		size_t count_valid_cluster = 0;
		size_t count_read = 0;
		while(test){
			printf("%zu - %zu : %zu reads // %s - %s \n",test->read->gene,test->read->mgene,test->pnbreads,chrNames[(int)(test->read->flags.chr)],chrNames[(int)(test->read->mchr)]);
			count_read+=test->pnbreads;
			count_valid_cluster ++;
			test=test->next;
		}

		size_t nb_total_valid = 0;
		MPI_Allreduce(&count_valid_cluster,&nb_total_valid,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
		if(!rank)
			printf("There are %zu VALID CLUSTERS\n",nb_total_valid);

		size_t nb_total_read = 0;
		MPI_Allreduce(&count_read,&nb_total_read,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
		if(!rank)
			printf("There are %zu READS\n",nb_total_read);


		//We set exons frame for each reads
		main_find_exon(&pClusters,count_valid_cluster,nbchr,gene_exon_path,chrNames,rank);

		//We do a some statistics for the result
		char* exon_valid_cluster = (char*)malloc(count_valid_cluster);
		test=pClusters;
		i=0;
		while(test){
			Read * r = test->read;
			size_t bool = 0;
			while(r){
				//We calculate how many reads have exon frame forwards and backwards which align perfectly
				//0->1,1->2,2->0
				if( (r->exon_coord >=0 && r->exon_mcoord >=0 && (r->exon_coord+1)%3==r->exon_mcoord)){
					bool++;
				}
				r=r->next;
			}
			exon_valid_cluster[i] = (bool*100)/test->pnbreads;
			i++;
			test=test->next;
		}
		tic=MPI_Wtime();

		//We write ours results
		fprintf(time_output_file, "%d ::::: *** WRITE ***\n", rank);
		writeCluster(rank, num_proc, pClusters, count_valid_cluster, chrNames, file_name, header, output_dir, exon_valid_cluster);
		fprintf(time_output_file, "%d (%.2lf)::::: *** FIN WRITE ***\n", rank, (double)(MPI_Wtime()-tic));

		free(exon_valid_cluster);
	}


	MPI_Barrier(MPI_COMM_WORLD);

	if(!rank){
		free(header);
	}

	free(readNumberByChr);

	for(i = 0; i < nbchr; i++)
		free(chrNames[i]);
	free(chrNames);

	fprintf(stderr,"Rank %d finished.\n", rank);
	MPI_Finalize();

	//Free file for time traces
	if(time_output_file)
		fclose(time_output_file);


	/////////////////////////TIME TO FREE///////////////////////////////

	while(pClusters){
		Cluster * cl_next = pClusters->next;
		Read * r = pClusters->read;
		while(r){
			Read * next = r->next;
			free(r);
			r=next;
		}
		free(pClusters);
		pClusters=cl_next;
	}
	free(reads);

	return 0;
}
