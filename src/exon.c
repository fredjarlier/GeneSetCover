/**
 * \file exon.c
 * \author Paul PANGANIBAN
 * \date August 26th 2015
 */

#include "exon.h"

void main_find_exon(Cluster ** pCluster,size_t nb_cluster,int nb_chr, char* gene_exon_path, char ** chrNames,int rank){

	int i = 0;
	char * exon_path = NULL;
	int ierr = 0;
	struct stat exon_file_stat;
	MPI_File mpi_exon_file;

	Cluster * c = *pCluster;
	Cluster_chain ** cl_by_chr =(Cluster_chain**)calloc(nb_chr-1,sizeof(Cluster_chain*));
	Cluster_chain ** first_cl_by_chr=(Cluster_chain**)calloc(nb_chr-1,sizeof(Cluster_chain*));
	Cluster_chain ** before_cl_by_chr=(Cluster_chain**)calloc(nb_chr-1,sizeof(Cluster_chain*));


	Cluster * cl_fake =(Cluster*)malloc(sizeof(Cluster));
	cl_fake->read=(Read*)malloc(sizeof(Read));
	cl_fake->read->gene=0;
	cl_fake->read->mgene=0;
	Read * r_fake =(Read*)malloc(sizeof(Read));
	r_fake->coord=0;
	r_fake->mcoord=0;

	//we reset our reads
	// coord and gene will be from forwards, then mcoord and mgene will be from backward
	while(c){
		Read * r = c->read;
		//printf("%zu-%zu :\n",c->read->gene,c->read->mgene);
		while(r){
			if(r->flags.is_mate == 0){
				size_t c_tmp = r->coord;
				r->coord =r->mcoord;
				r->mcoord=c_tmp;

				unsigned char chr_tmp = r->mchr;
				r->mchr=r->flags.chr;
				r->flags.chr=chr_tmp;
			}

			if( r->flags.replace_gene_with_mgene == r->flags.is_mate){
				size_t tmp=r->gene;
				r->gene = r->mgene;
				r->mgene=tmp;
			}
			//printf("\t%zu-%zu: %zu-%zu // %s & %s // IM : %u // RG :%u\n ",r->gene,r->mgene,r->coord,r->mcoord,chrNames[r->flags.chr],chrNames[r->mchr],r->flags.is_mate,r->flags.replace_gene_with_mgene);
			r=r->next;
		}
		//we trie reads by coord to have forward next of backward
		r_fake->next=c->read;
		mergeSort(r_fake,c->pnbreads);
		c->read=r_fake->next;
		c=c->next;
	}

	free(r_fake);
	c=(*pCluster);


	//first we sort by gene
	if(c && c->next){
		cl_fake->next=c;
		mergeSortCluster_by_Gene(cl_fake,nb_cluster);
		c=cl_fake->next;
	}
	(*pCluster)=c;


	//we group ours clusters by chr in struct cluster_chain
	while(c){
		Cluster_chain * clc_tmp =(Cluster_chain*)malloc(sizeof(Cluster_chain));

		unsigned char chr_forward = c->read->flags.chr;
		//we want to open forwards'chr

		clc_tmp->cl=c;
		clc_tmp->next = NULL;
		if(before_cl_by_chr[chr_forward] == NULL){
			before_cl_by_chr[chr_forward]= clc_tmp;
			first_cl_by_chr[chr_forward]=  clc_tmp;
		}
		else{
			before_cl_by_chr[chr_forward]->next=clc_tmp;
			before_cl_by_chr[chr_forward] = before_cl_by_chr[chr_forward]->next;
		}
		c=c->next;
	}

	//we go back to first
	for(i=0;i<nb_chr-1;i++){
		cl_by_chr[i]=first_cl_by_chr[i];
	}


	//we open each gene_exon files
	for(i=0;i<nb_chr-1;i++){
		exon_path=(char*)calloc(255,sizeof(char));
		strcpy(exon_path,gene_exon_path);
		strcat(exon_path,chrNames[i]);
		strcat(exon_path,"_exon.gtf");

		ierr = MPI_File_open(MPI_COMM_WORLD, exon_path,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_exon_file);
		if (ierr){
			if(!rank)
			DEBUG("Failed to open file in process 0 %s\n",exon_path);
		}

		//we take file size
		else{
			int file_size = stat(exon_path, &exon_file_stat);
			if (file_size == -1){
				fprintf(stderr,"Failed to find file size\n");
				MPI_Abort(MPI_COMM_WORLD, 304);
				exit(0);
			}
			file_size = (long long)exon_file_stat.st_size;

			if(file_size){
				char * buffer =(char*)calloc(file_size+1,sizeof(char));
				MPI_Status chr_status;
				MPI_File_read_all(mpi_exon_file,buffer,file_size,MPI_CHAR,&chr_status);
				get_exonForward_for_each_cluster(cl_by_chr[i],buffer);
				free(buffer);
			}

		}
		MPI_File_close(&mpi_exon_file);
		free(exon_path);
	}


	//we have forwards exon for each reads, now we do the same for backwards
	for(i=0;i<nb_chr-1;i++){
		while(cl_by_chr[i] && cl_by_chr[i]->next){
			Cluster_chain * tmp_next = cl_by_chr[i]->next->next;
			free(cl_by_chr[i]->next);
			cl_by_chr[i]->next=tmp_next;
		}
		if(i!=0 && cl_by_chr[i])
			free(cl_by_chr[i]);
	}
	free(cl_by_chr);
	free(before_cl_by_chr);
	free(first_cl_by_chr);

	cl_by_chr =(Cluster_chain**)calloc(nb_chr-1,sizeof(Cluster_chain*));
	first_cl_by_chr=(Cluster_chain**)calloc(nb_chr-1,sizeof(Cluster_chain*));
	before_cl_by_chr=(Cluster_chain**)calloc(nb_chr-1,sizeof(Cluster_chain*));
	c=(*pCluster);

	//first we sort by gene
	if(c && c->next){
		cl_fake->next=c;
		mergeSortCluster_by_Mgene(cl_fake,nb_cluster);
		c=cl_fake->next;
	}
	(*pCluster)=c;

	//we group ours clusters by chr in struct cluster_chain
	while(c){
		Cluster_chain * clc_tmp =(Cluster_chain*)malloc(sizeof(Cluster_chain));
		//if(rank==46)
		//printf("\t%zu-%zu\n",c->read->gene,c->read->mgene);
		unsigned char chr_backward = c->read->mchr;
		//we want to open forwards'chr


		clc_tmp->cl=c;
		clc_tmp->next = NULL;
		if(before_cl_by_chr[chr_backward] == NULL){
			before_cl_by_chr[chr_backward]= clc_tmp;
			first_cl_by_chr[chr_backward]=  clc_tmp;
		}
		else{
			before_cl_by_chr[chr_backward]->next=clc_tmp;
			before_cl_by_chr[chr_backward] = before_cl_by_chr[chr_backward]->next;
		}
		c=c->next;
	}


	//we go back to first
	for(i=0;i<nb_chr-1;i++){
		cl_by_chr[i]=first_cl_by_chr[i];
	}


	//we open each gene_exon files
	for(i=0;i<nb_chr-1;i++){
		exon_path=(char*)calloc(255,sizeof(char));
		strcpy(exon_path,gene_exon_path);
		strcat(exon_path,chrNames[i]);
		strcat(exon_path,"_exon.gtf");

		ierr = MPI_File_open(MPI_COMM_WORLD, exon_path,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_exon_file);
		if (ierr){
			if(!rank)
			DEBUG("Failed to open file in process 0 %s\n",exon_path);
		}

		else{
			//we take file size
			int file_size = stat(exon_path, &exon_file_stat);
			if (file_size == -1){
				fprintf(stderr,"Failed to find file size\n");
				MPI_Abort(MPI_COMM_WORLD, 304);
				exit(0);
			}
			file_size = (long long)exon_file_stat.st_size;

			if(file_size){
				char * buffer =(char*)calloc(file_size+1,sizeof(char));
				MPI_Status chr_status;
				MPI_File_read_all(mpi_exon_file,buffer,file_size,MPI_CHAR,&chr_status);
				get_exonBackward_for_each_cluster(cl_by_chr[i],buffer);
				free(buffer);
			}
		}

		MPI_File_close(&mpi_exon_file);
		free(exon_path);
	}

	for(i=0;i<nb_chr-1;i++){
		while(cl_by_chr[i] && cl_by_chr[i]->next){
			Cluster_chain * tmp_next = cl_by_chr[i]->next->next;
			free(cl_by_chr[i]->next);
			cl_by_chr[i]->next=tmp_next;
		}
		if(i!=0 && cl_by_chr[i])
			free(cl_by_chr[i]);
	}

	free(cl_fake);
	free(cl_by_chr);
	free(before_cl_by_chr);
	free(first_cl_by_chr);
}


void get_exonForward_for_each_cluster(Cluster_chain * cl_by_chr, char * buffer_file){
	Cluster_chain * clc = cl_by_chr;
	int i = 0 ;
	size_t size_buffer=strlen(buffer_file);
	int size_line =0;
	size_t id_curr = 0;
	size_t gene_curr=0;
	while(clc && i<size_buffer){
		//first we find the line for our gene
		size_line =0;

		gene_curr=clc->cl->read->gene;

		//we find the id
		while(buffer_file[i+size_line] != '\t' && i+size_line<size_buffer){
			size_line++;
		}
		buffer_file[i+size_line]=0;
		id_curr=strtoul(buffer_file+i,NULL,10);

		size_line++;



		//if its the good line we check where is our transcript
		if(id_curr == gene_curr){

			int nb_exon=0;
			int p_nb_exon = size_line;

			//we prepare nb_exon
			char * gene_forward = buffer_file+i+size_line;
			while(buffer_file[i+size_line] != '\t' && i+size_line<size_buffer){
				size_line++;
			}
			buffer_file[i+size_line] = 0;
			size_line++;

			//we prepare nb_exon
			while(buffer_file[i+size_line] != '\t' && i+size_line<size_buffer){
				size_line++;
			}
			buffer_file[i+size_line] = 0;
			nb_exon=strtoul(buffer_file+i+p_nb_exon,NULL,10);

			size_line++;

			char * start_exon = buffer_file+i+size_line;

			while(buffer_file[i+size_line]!='\t'){
				if(buffer_file[i+size_line]==','){
					buffer_file[i+size_line]=0;
				}
				size_line++;
			}
			size_line++;


			char * end_exon = buffer_file+i+size_line;
			while(buffer_file[i+size_line]!='\t'){
				if(buffer_file[i+size_line]==','){
					buffer_file[i+size_line]=0;
				}
				size_line++;
			}
			size_line++;
			char * exon_frame = buffer_file+i+size_line;

			while(buffer_file[i+size_line]!='\n' && (i+size_line)<size_buffer){
				if(buffer_file[i+size_line]==','){
					buffer_file[i+size_line]=0;
				}
				size_line++;
			}
			size_line++;

			while(clc && id_curr == gene_curr){
				get_exon_forward_for_each_reads(clc->cl,start_exon,end_exon,exon_frame,nb_exon);
				clc->cl->gene_name=(char*)calloc(strlen(gene_forward)+1,sizeof(char));
				strcpy(clc->cl->gene_name,gene_forward);

				clc=clc->next;
				if(clc && clc->cl){
					gene_curr=clc->cl->read->gene;
				}
			}
			i+=size_line;//start of new line
		}


		//else we look next line
		else{
			while(buffer_file[i+size_line]!='\n' && (i+size_line)<size_buffer){
				size_line++;
			}
			size_line++;
			i+=size_line;
		}

	}//while clc and i

	if(clc){
		printf("ID GENE AND GENE NOT FIND %zu : %zu-%zu / %u-%u !!!!!!!!!!! \n",gene_curr,clc->cl->read->gene,clc->cl->read->mgene,clc->cl->read->flags.chr,clc->cl->read->mchr);
		MPI_Abort(MPI_COMM_WORLD,gene_curr);
	}
}


void get_exon_forward_for_each_reads(Cluster*clust,char* start_exon, char* end_exon, char* exon_frame,int nb_exon_total){

	int nb_exon_skiped = 0;
	int j_start_exon=0;
	int j_end_exon=0;
	int j_exon_frame=0;
	Read * reads = clust->read;
	Read* r_curr = reads;
	size2_t frame;
	size_t coord_curr=0;

	//first we checked our first coord
	coord_curr = r_curr->coord;

	while(r_curr && nb_exon_skiped<nb_exon_total){
		//we place our end of frame
		frame.b=strtoul(end_exon+j_end_exon,NULL,10);
		while(frame.b < coord_curr && nb_exon_skiped<nb_exon_total){
			while(end_exon[j_end_exon]!=0){
				j_end_exon++;
			}
			j_end_exon++;
			nb_exon_skiped++;
			if(nb_exon_skiped<nb_exon_total){
				frame.b=strtoul(end_exon+j_end_exon,NULL,10);
			}
		}

		//2 cases : out of bounds or find a frame

		//Case 1 : out of bounds
		//all reads have exon set to 3(not found)


		if(nb_exon_skiped>=nb_exon_total){
			while(r_curr){
				r_curr->exon_coord=(-2);
				r_curr=r_curr->next;
			}

		}

		//Case 2 : we look if reads are in one frame or not
		else{
			int i =0;
			//we go to the exon start of our frame found
			for(i=0;i<nb_exon_skiped;i++){
				while(start_exon[j_start_exon]!=0){
					j_start_exon++;
				}
				j_start_exon++;
			}
			frame.a=strtoul(start_exon+j_start_exon,NULL,10);

			//we do it for all reads that are in this frame
			while(r_curr && coord_curr<frame.b){
				int exon = (-2);
				//if coord is not in the frame
				if(frame.a>coord_curr+50){
					r_curr->exon_coord=exon;
				}
				//this read match on this exon
				else{
					for(i=0;i<nb_exon_skiped;i++){
						while(exon_frame[j_exon_frame]!=0){
							j_exon_frame++;
						}
						j_exon_frame++;
					}
					exon = atoi(exon_frame);
					r_curr->exon_coord=exon;
				}
				r_curr=r_curr->next;
				if(r_curr){
					coord_curr=r_curr->coord;
				}
			}
		}

		//if there are still reads , we look in an other frame
		if(r_curr){
			while(end_exon[j_end_exon]!=0 && nb_exon_skiped<nb_exon_total){
				j_end_exon++;
			}
			j_end_exon++;
			nb_exon_skiped++;
		}
	}

}


void get_exonBackward_for_each_cluster(Cluster_chain * cl_by_chr, char * buffer_file){
	Cluster_chain * clc = cl_by_chr;
	int i = 0 ;
	size_t size_buffer=strlen(buffer_file);
	int size_line =0;
	size_t id_curr = 0;
	size_t gene_curr=0;
	while(clc && i<size_buffer){
		//first we find the line for our gene
		size_line =0;

		gene_curr=clc->cl->read->mgene;

		//we find the id
		while(buffer_file[i+size_line] != '\t' && i+size_line<size_buffer){
			size_line++;
		}
		buffer_file[i+size_line]=0;
		id_curr=strtoul(buffer_file+i,NULL,10);

		size_line++;



		//if its the good line we check where is our transcript
		if(id_curr == gene_curr){


			char * gene_backward = buffer_file+i+size_line;
			while(buffer_file[i+size_line] != '\t' && i+size_line<size_buffer){
				size_line++;
			}
			buffer_file[i+size_line] = 0;
			size_line++;


			int nb_exon=0;
			int p_nb_exon = size_line;

			//we prepare nb_exon
			while(buffer_file[i+size_line] != '\t' && i+size_line<size_buffer){
				size_line++;
			}
			buffer_file[i+size_line] = 0;
			nb_exon=strtoul(buffer_file+i+p_nb_exon,NULL,10);

			size_line++;

			char * start_exon = buffer_file+i+size_line;

			while(buffer_file[i+size_line]!='\t'){
				if(buffer_file[i+size_line]==','){
					buffer_file[i+size_line]=0;
				}
				size_line++;
			}
			size_line++;


			char * end_exon = buffer_file+i+size_line;
			while(buffer_file[i+size_line]!='\t'){
				if(buffer_file[i+size_line]==','){
					buffer_file[i+size_line]=0;
				}
				size_line++;
			}
			size_line++;

			char * exon_frame = buffer_file+i+size_line;

			while(buffer_file[i+size_line]!='\n' && (i+size_line)<size_buffer){
				if(buffer_file[i+size_line]==','){
					buffer_file[i+size_line]=0;
				}
				size_line++;
			}
			size_line++;

			while(clc && id_curr == gene_curr){
				get_exon_backward_for_each_reads(clc->cl,start_exon,end_exon,exon_frame,nb_exon);
				clc->cl->mgene_name=(char*)calloc(strlen(gene_backward)+1,sizeof(char));
				strcpy(clc->cl->mgene_name,gene_backward);
				clc=clc->next;
				if(clc && clc->cl){
					gene_curr=clc->cl->read->mgene;
				}
			}
			i+=size_line;//start of new line
		}


		//else we look next line
		else{
			while(buffer_file[i+size_line]!='\n' && (i+size_line)<size_buffer){
				size_line++;
			}
			size_line++;
			i+=size_line;
		}

	}//while clc and i

	if(clc){
		printf("ID GENE AND MGENE NOT FIND %zu : %zu-%zu / %u-%u !!!!!!!!!!! \n",gene_curr,clc->cl->read->gene,clc->cl->read->mgene,clc->cl->read->flags.chr,clc->cl->read->mchr);
		MPI_Abort(MPI_COMM_WORLD,gene_curr);
	}
}


void get_exon_backward_for_each_reads(Cluster*clust,char* start_exon, char* end_exon, char* exon_frame,int nb_exon_total){

	int nb_exon_skiped = 0;
	int j_start_exon=0;
	int j_end_exon=0;
	int j_exon_frame=0;
	Read * reads = clust->read;
	Read* r_curr = reads;
	size2_t frame;
	size_t coord_curr=0;

	//first we checked our first coord
	coord_curr = r_curr->mcoord;

	while(r_curr && nb_exon_skiped<nb_exon_total){
		//we place our end of frame
		frame.b=strtoul(end_exon+j_end_exon,NULL,10);
		while(frame.b < coord_curr && nb_exon_skiped<nb_exon_total){
			while(end_exon[j_end_exon]!=0){
				j_end_exon++;
			}
			j_end_exon++;
			nb_exon_skiped++;
			if(nb_exon_skiped<nb_exon_total){
				frame.b=strtoul(end_exon+j_end_exon,NULL,10);
			}
		}

		//2 cases : out of bounds or find a frame

		//Case 1 : out of bounds
		//all reads have exon set to 3(not found)

		//printf("\tnb_exon_skiped : %d\n",nb_exon_skiped);

		if(nb_exon_skiped>=nb_exon_total){
			while(r_curr){
				r_curr->exon_mcoord=(-2);
				r_curr=r_curr->next;
			}
		}

		//Case 2 : we look if reads are in one frame or not
		else{
			int i =0;
			//we go to the exon start of our frame found
			for(i=0;i<nb_exon_skiped;i++){
				while(start_exon[j_start_exon]!=0){
					j_start_exon++;
				}
				j_start_exon++;
			}
			frame.a=strtoul(start_exon+j_start_exon,NULL,10);

			//we do it for all reads that are in this frame
			while(r_curr && coord_curr<frame.b){
				int exon = (-2);
				//if coord is not in the frame
				if(frame.a>coord_curr+50){
					r_curr->exon_mcoord=exon;
				}
				//this read match on this exon
				else{
					for(i=0;i<nb_exon_skiped;i++){
						while(exon_frame[j_exon_frame]!=0){
							j_exon_frame++;
						}
						j_exon_frame++;
					}
					exon = atoi(exon_frame);
					r_curr->exon_mcoord=exon;
				}
				r_curr=r_curr->next;
				if(r_curr){
					coord_curr=r_curr->mcoord;
				}
			}
		}

		//if there are still reads , we look in an other frame
		if(r_curr){
			while(end_exon[j_end_exon]!=0 && nb_exon_skiped<nb_exon_total){
				j_end_exon++;
			}
			j_end_exon++;
			nb_exon_skiped++;
		}
	}

}
