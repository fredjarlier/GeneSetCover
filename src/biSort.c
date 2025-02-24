/**
 * \file biSort.c
 * \author Thomas MAGALHAES, Paul PANGANIBAN, Nizar HDADECH
 * \date August 26th 2015
 */

#include "biSort.h"
#include "mergeSort.h"

size_t biSort(Read ***reads,Cluster** pCluster,size_t *preadNumberByChr,int rank,size_t * nbCluster)
{	
	Cluster* clusters=NULL;
	Cluster * first = NULL;
	size_t gene , mgene;
	size_t preadNumber=0;

	Read* anchor = NULL;//first read of the new linked list
	Read* before = NULL;//
	int i=0,j=1;
	int first_chr = -1;
	Read* before_read_chr = NULL;

	//MERGE ALL READS BY CHR IN ONE LIST OF READS
	Read * r_curr = NULL;
	for(i=0;i<25;i++){
		r_curr=(*reads)[i];

		for(j=0;j<preadNumberByChr[i]-1 && preadNumberByChr[i]>0;j++){
			if(r_curr && r_curr->gene > r_curr->mgene){
				size_t gene = r_curr->gene;
				r_curr->gene=r_curr->mgene;
				r_curr->mgene=gene;
				r_curr->flags.replace_gene_with_mgene=1; //use to know if gene is from coor or mcoord
			}
			r_curr = r_curr->next;
		}


		if(first_chr <0 && !before_read_chr && preadNumberByChr[i]>0){
			first_chr = i ;
			before_read_chr = r_curr;
		}

		else if(before_read_chr && preadNumberByChr[i]>0 && i<24){
			before_read_chr ->next = (*reads)[i];
			before_read_chr = r_curr;
		}
	}


	//if we have reads
	if(first_chr>=0){
		r_curr=(*reads)[first_chr]; //we go to our first read
		anchor = r_curr;

		//first we find our before
		while(before == NULL && r_curr !=NULL){
			//if it s discordant we keep it
			if(r_curr->gene != 0 && r_curr->gene!=r_curr->mgene && r_curr->mgene !=0){
				preadNumber=1;
				before = r_curr;
				anchor = before;
				r_curr = r_curr->next;
			}

			//we free it
			else{
				anchor = r_curr->next;
				free(r_curr);
				r_curr=anchor;
				if(r_curr != NULL){
					r_curr=r_curr->next;
				}
			}
		}

		//CONDITION 1 on all reads
		while(r_curr){
			//We keep it if it's discordant
			if(r_curr->gene != 0 && r_curr->gene!=r_curr->mgene && r_curr->mgene !=0 ){
				preadNumber++;
				before = r_curr;
				r_curr=r_curr->next;
			}

			else{
				//we delete it
				before->next = r_curr->next;
				free(r_curr);
				r_curr=before->next;
			}
		}


		//ON REVIENT AU PREMIER READ
		r_curr = anchor;

		size_t count =0;
		while(r_curr){
			count++;
			r_curr=r_curr->next;
		}

		preadNumber=count;
		r_curr = anchor;

		if(r_curr){
			anchor=(Read*)malloc(sizeof(Read));
			anchor->coord=1;
			anchor->mcoord=2;
			anchor->gene=1;
			anchor->mgene=2;
			anchor->flags.chr=0;
			anchor->flags.replace_gene_with_mgene=0;
			anchor->flags.is_mate=0;
			anchor->flags.left=0;
			anchor->link=NULL;
			anchor->offset=0;
			anchor->offset_source_file=0;
			anchor->quality=0;

			anchor->next = r_curr;
			r_curr=anchor;
			preadNumber++;

			//WE DO THE BI SORT
			if(r_curr->next && r_curr->next->next)
				mergeSortGene(r_curr,preadNumber);

			r_curr=r_curr->next;
			free(anchor);
			preadNumber--;
			anchor = r_curr;


			//we init the first cluster
			clusters =(Cluster*)malloc(sizeof(Cluster));
			first = clusters;
			clusters->read=r_curr;
			clusters->pnbreads=1;
			clusters->next=NULL;
			clusters->process_rank=rank;
			gene = r_curr->gene;
			mgene = r_curr->mgene;
			(*nbCluster)++;

			Read * r_before = r_curr;
			r_curr = r_curr->next;
			for(i=0;i<(preadNumber-1);i++)
			{
				//if it's the same cluster we add it
				if(r_curr->gene == gene && r_curr->mgene==mgene){
					clusters->pnbreads++;
					r_before=r_curr;
					r_curr=r_curr->next;
				}
				//if it s a new cluster
				else{
					r_before->next=NULL;
					//we create it
					clusters->next = (Cluster*)malloc(sizeof(Cluster));
					clusters = clusters->next;
					clusters->read=r_curr;
					clusters->pnbreads=1;
					clusters->next=NULL;
					clusters->process_rank=rank;
					gene = r_curr->gene;
					mgene = r_curr->mgene;
					(*nbCluster)++;

					r_before = r_curr;
					r_curr=r_curr->next;
				}

			}
		}
	}

	*pCluster=first;
	(**reads)=anchor;

	return preadNumber;
}


void getGeneFamily(char * buffer_family, Cluster * pCluster)
{
	int i = 0;
	int j =0;
	int size_line=0;

	Cluster * cl = pCluster;

	size_t gene = cl->read->gene;
	size_t size_buff =strlen(buffer_family);
	while( i<size_buff && cl){

		size_line=0;
		j=0;

		//we find the size of current line and copy it into current_line
		while(buffer_family[i+size_line] != '\n' && buffer_family[i+size_line] != 0 ){
			size_line++;
		}
		buffer_family[i+size_line]=0;

		//now we take the id
		j=0;
		while(buffer_family[i+j]!='\t'){
			j++;
		}
		buffer_family[i+j]=0;
		size_t curr_gene = strtoul(buffer_family+i,NULL,10);
		j++;
		//if its the good gene , we look the family id
		if(curr_gene==gene){
			size_t family = strtoul(buffer_family+i+j,NULL,10);
			//we add this family id to the cl->next with the same gene
			//thanks to bisort
			while(cl && gene == curr_gene){
				cl->family_gene= family;
				cl = cl->next;
				if(cl){
					gene=cl->read->gene;
				}
				else{
					gene = 0;
				}
			}
		}

		//we move the seek
		i+=(size_line+1);
	}
}


size_t getMGeneFamily(char * buffer_family, Cluster ** pCluster)
{
	int i = 0;
	int j =0;
	int size_line=0;

	size_t nb_same_family = 0;

	//skip first line

	Cluster * cl = *pCluster;
	Cluster * cl_before = NULL;

	size_t mgene = cl->read->mgene;
	size_t size_buff =strlen(buffer_family);
	while( i<size_buff && cl){

		size_line=0;
		j=0;

		//we find the size of current line and copy it into current_line
		while(buffer_family[i+size_line] != '\n' && buffer_family[i+size_line] != 0 ){
			size_line++;
		}
		buffer_family[i+size_line]=0;

		//now we take the id
		j=0;
		while(buffer_family[i+j]!='\t'){
			j++;
		}
		buffer_family[i+j]=0;
		size_t curr_mgene = strtoul(buffer_family+i,NULL,10);
		j++;
		//if its the good gene , we look the family id
		if(curr_mgene==mgene){
			size_t family = strtoul(buffer_family+i+j,NULL,10);
			//we add this family id to the cl->next with the same gene
			//thanks to bisort
			while(cl && mgene == curr_mgene){
				cl->family_mgene= family;

				if(cl->family_gene == cl->family_mgene){
					nb_same_family++;

					//free reads
					Read * r_next=cl->read->next;
					while(r_next){
						cl->read->next=r_next->next;
						free(r_next);
						r_next=cl->read->next;
					}
					free(cl->read);

					//free cluster
					if(cl_before){
						cl_before->next=cl->next;
						free(cl);
						cl=cl_before->next;
					}
					else{
						(*pCluster) = cl->next;
						free(cl);
						cl=(*pCluster);
					}
				}

				//if fusion between different familys
				else{
					cl_before=cl;
					cl = cl->next;
				}

				if(cl){
					mgene=cl->read->mgene;
				}
				else{
					mgene = 0;
				}
			}
		}
		//we move the seek
		i+=(size_line+1);
	}
	return(nb_same_family);
}


