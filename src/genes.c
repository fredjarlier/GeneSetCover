/**
 * \file genes.c
 * \author Thomas MAGALHAES, Paul PANGANIBAN
 * \date August 26th 2015
 */

#include "genes.h"
#include "assert.h"


size2_t getGeneAsNumber(char* buffer_file, Read* read){

	size_t file_length = strlen(buffer_file);
	size_t i = 0;
	size2_t found ;
	Read * r = read->next;
	Read * r_before = read;

	size_t result = 0;
	int size_line = 0;
	size2_t size;


	found.a=0;
	found.b=0;
	int bool_next_gene = 0;

	while(i<file_length && r){

		//find the size of the line
		size_line=0;
		while(buffer_file[i+size_line] != '\n' && buffer_file[i+size_line] != 0 ){
			size_line++;
		}

		buffer_file[i+size_line] = 0;

		//return the gene value and the bounds of this gene
		result = getGeneLineBounds(buffer_file+i, &size);
		bool_next_gene =0;
		while(!bool_next_gene && r!=NULL){
			if(size.a == -1 || size.b == -1)
				{
					printf("EROR IN GENE FILE\n");
					i+=(size_line+1);
					bool_next_gene = 1;
					continue;
				}


				//if coord is in the position frame
				//we set gene with the gene id in that line
				if((r->coord+50) >= (size.a))
				{
					if(r->coord <= size.b)
					{
						found.a++;
						r->gene = result;

						r_before = r;
						r = r->next;
					}
					else{
						i+=(size_line+1);
						bool_next_gene=1;
					}
				}


				//if its coord is out of bounds we delete it
				else
				{
					r_before ->next = r->next;
					free(r);
					r=r_before->next;
					found.b++;
				}
		}//while bool_next_gne

		//if there is an error in gene file
			//free(current_line);
	}
	//if there are still reads , we delete them because they are out of bounds
	while(r){
		r_before->next=r->next;
		free(r);
		r=r_before->next;
		found.b++;
	}
	return found;
}


size2_t getMGeneAsNumber(char* buffer_file, Read_chain* backward){

	size_t file_length = strlen(buffer_file);
	size_t i = 0;
	size2_t found ;
	Read_chain * r = backward;

	size_t result = 0;
	int size_line = 0;
	size2_t size;


	found.a=0;
	found.b=0;

	int bool_next_gene = 0;

	while(i<file_length && r->reads){

		//find the size of the line
		size_line=0;
		while(buffer_file[i+size_line] != '\n' && buffer_file[i+size_line] != 0 ){
			size_line++;
		}
		buffer_file[i+size_line] = 0;

		//return the gene value and the bounds of this gene
		result = getGeneLineBounds(buffer_file+i, &size);
		bool_next_gene = 0;

		while(!bool_next_gene && r->reads){
			//if there is an error in gene file
				if(size.a == -1 || size.b == -1)
				{
					printf("EROR IN GENE FILE\n");
					i+=(size_line+1);
					bool_next_gene=1;
					continue;
				}


				//if mcoord is in the position frame
				//we set gene with the mgene id in that line
				if((r->reads->mcoord+50) >= (size.a))
				{
					if(r->reads->mcoord <= size.b)
					{
						r->reads->mgene = result;
						r = r->next;
						found.a++;
					}
					//else we go to next line
					else{
						i+=(size_line+1);
						bool_next_gene=1;
					}
				}
				//if its coord is out of bounds we check next backward
				else
				{
					if(!(r->reads->flags.is_mate))
						r->reads->mgene = 0;
					else
						r->reads->gene = 0;

					r = r->next;
					found.b++;
				}
		}//while bool_next_gene
	}

	//Reads out of bounds
	while(r->reads){
		if(!(r->reads->flags.is_mate))
			r->reads->mgene = 0;
		else
			r->reads->gene = 0;
		r=r->next;
		found.b++;
	}
	return found;
}


size_t getGeneLineBounds(char* line, size2_t* size)
{
	char* buf = line;
	int count, a = 0, b = 0;
	int i, j;
	size_t value = 0;
	value = -1;
	size->a = 0;
	size->b = 0;

	count = strlen(buf);

	//4e ',' = start pos
	j = 0;
	for(i = 0; i<count; i++)
	{
		//if(buf[i] == ',')
		if(buf[i] == '\t')
		{
			++j;

			if(j == 2)
			{
				buf[i] = '\0';
				value = ++i;
				break;
			}
			else if(j == 0)
			{
				buf[i] = '\0';
				a = ++i;
			}
			else if(j == 1)
			{
				buf[i] = '\0';
				b = ++i;

			}//if - else if
		}//if
	}//for

	char *cb =  buf + a;
	char *cb1 =  buf + b;
	char *cb2 =  buf + value ;
	size->a = strtoul( cb, &cb, 10);
	size->b = strtoul( cb1, &cb1, 10);
	value = strtoul( cb2, &cb2, 10);

	if(size->a > size->b)
	{
		fprintf(stderr, "Error with genes bounds\n");
		size->a = -1;
		size->b = -1;
	}


	return value;
}





Read_chain** getBackward(Read** reads, int nbchr, size_t* count)
{
	int i;
	Read_chain** backward;
	if(!count)
		count = (size_t*)calloc(nbchr,sizeof(size_t));

	Read_chain* tmp_back;
	Read* tmp_read = NULL;
	unsigned char tmp_chr = 0;
	backward=(Read_chain**)malloc(sizeof(Read_chain*)*nbchr);
	assert(backward != NULL);
	//Init
	for(i = 0; i < nbchr-1; i++){
		count[i] = 0;

		backward[i] = (Read_chain*)malloc(sizeof(Read_chain));
		backward[i]->reads = NULL;
		backward[i]->next = NULL;
		assert(backward[i] != NULL);
	}

	//Set
	for(i = 0; i < nbchr-1; i++){

		if(reads[i])
			tmp_read = reads[i];

		while(tmp_read != NULL)
		{
			tmp_chr = 0;
			//assert((tmp_read->flags.chr >= 0) && (tmp_read->flags.chr <= 25));
			tmp_chr = tmp_read->mchr;
			
			if(tmp_chr > (unsigned char)(nbchr - 1))
				tmp_chr = 0;
			
			//tmp_back = backward[tmp_chr];

			tmp_back = (Read_chain*)malloc(sizeof(Read_chain));
			assert(tmp_back != NULL);
			tmp_back->reads = NULL;
			tmp_back->next = NULL;

			tmp_back->next = backward[tmp_chr];
			backward[tmp_chr] = tmp_back;
			tmp_back->reads = tmp_read;

			tmp_read = tmp_read->next;
			count[tmp_chr]++;
		}
	}
	assert(backward != NULL);
	return backward;
}


void freeBackward(Read_chain** chain, int nbchr)
{
	int i=0;
	Read_chain* tmp=NULL;

	for(i = 0; i < nbchr-1; i++)
	{
		//freeSubBackward(chain[i]);
		//tmp = chain[i]->next;
		while(chain[i] != NULL)
		{
			tmp = chain[i]->next;
			free(chain[i]);
			chain[i]=tmp;
		}
	}
	free(chain);
}





