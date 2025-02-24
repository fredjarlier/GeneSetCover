/**
 * \file parser.c
 * \author Nicolas FEDY, Leonor SIROTTI, Paul PANGANIBAN, Thomas MAGALHAES
 * \date August 26th 2015
 */

#include "parser.h"

unsigned int find_header(char *localData, int rank, size_t *unmappedSize, int *pnbchr, char **pheader, char ***pchrNames){
	char *currentCarac;
	char *header;
	char **chrNames = NULL;
	char currentLine[MAX_LINE_SIZE];
	unsigned int i, nbchr = 0;
	int next;
	size_t headerSize, lineSize;
	size_t j;

	//we take the first line
	next = tokenizer(localData,'\n', currentLine);

	headerSize = 0; //offset size of header to replace offset of proc 0
	lineSize = 0;//use to update headerSize

	chrNames = (char**)malloc(sizeof(char*));

	//RANK 0 WILL CALCULATE SIZE HEADER AND FIND CHRNAME
	if(rank == 0){
		i = 0;
		j = 0;
		header = (char*)malloc(sizeof(char));
		*header = 0;

		//we look headers lines
		//we add header string into variable 'header' and add find chrNames
		while (next && currentLine[0] == '@'){

			currentCarac = currentLine;

			lineSize = strlen(currentLine) + 1;
			headerSize += lineSize;

			header = realloc(header, (headerSize+1)*sizeof(char));

			//add currentline into header
			while(*currentCarac){
				header[j++] = *currentCarac++;
			}

			header[j++] = '\n';
			header[j] = 0;
			*pheader = header;

			//find chromosome name
			if(currentLine[1] == 'S' && currentLine[2] == 'Q'){

				currentCarac = strtok(currentLine, "\t");
				currentCarac = strtok(NULL, "\t");
				currentCarac += 3;

				chrNames[i] = (char*)malloc((strlen(currentCarac)+1)*sizeof(char));

				strcpy(chrNames[i], currentCarac);

				nbchr++;
				chrNames = (char**)realloc(chrNames, (nbchr+1)*sizeof(char*));

				i++;
			}

			next = tokenizer(NULL, '\n', currentLine);
		}

		//we add UNMAPPED to chrNames
		chrNames[i] = (char*)malloc((strlen(UNMAPPED)+1)*sizeof(char));
		strcpy(chrNames[i], UNMAPPED);
		nbchr++;
	}

	/*************************
	 * GIVE THE HEADER SIZE TO ALL
	 *************************/
	if(rank==0){
		printf("HEADERSIZE : %zu\n",headerSize);
		printf("%s\n",header);
		printf("NB CHR : %d\n",nbchr);
	}
	MPI_Bcast(&headerSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	*unmappedSize = headerSize;


	/**********************
	 * GIVE CHR NAMES TO ALL
	 *************************/

	//1 - every process will know nbchr
	MPI_Bcast(&nbchr, 1, MPI_INT, 0,MPI_COMM_WORLD);
	*pnbchr = nbchr;

	if(nbchr!=0){
		size_t size_chrName[nbchr];//list of size of chromosome name
		size_t size_all_chrname = 0;
		for(i=0;i<nbchr;i++){
			size_chrName[i]=0;
		}

		//2 - rank 0 will set the size of each chrname
		if(rank ==0){
			for(i=0;i<nbchr;i++){
				size_chrName[i]=strlen(chrNames[i]);
				size_all_chrname+=size_chrName[i];
			}
			size_all_chrname++;
		}

		//3 - Other process will prepare buff for each chromosome name
		MPI_Bcast(size_chrName,nbchr,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
		if(rank){
			chrNames = (char**)malloc((nbchr)*sizeof(char*));

			for(i=0;i<nbchr;i++){
				chrNames[i]=(char*)malloc(size_chrName[i]*sizeof(char)+1);
				*(chrNames[i]) = 0;
				size_all_chrname += size_chrName[i];
			}
			size_all_chrname++;
		}

		//4 - We will send all read names in one buffer
		char * buff_chrNames =(char*)malloc(size_all_chrname*sizeof(char));
		for(i=0;i<size_all_chrname;i++){
			buff_chrNames[i]=0;
		}

		if(rank==0){
			for(i=0;i<nbchr;i++){
				strcat(buff_chrNames,chrNames[i]);
			}
		}

		//5 - Other rank will add chr names
		MPI_Bcast(buff_chrNames,(int)size_all_chrname,MPI_CHAR,0,MPI_COMM_WORLD);
		if(rank){
			int offset_chrname = 0;
			for(i=0;i<nbchr;i++){
				int n = 0;
				for(n=0;n<size_chrName[i];n++){
					chrNames[i][n]=buff_chrNames[offset_chrname+n];
				}
				chrNames[i][size_chrName[i]]='\0';
				offset_chrname+=size_chrName[i];
			}
		}

		free(buff_chrNames);
		*pchrNames = chrNames;
	}
	return(headerSize);
}

size_t * init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank){
	size_t * goff =(size_t*)calloc((size_t)(numproc+1), sizeof(size_t));
	char * current_line = NULL;
	MPI_Status status;
	int i = 0;
	int j = 0;

	size_t lsize = fsize/numproc;
	goff[0]=headerSize;
	for(i=1;i<numproc;i++){
		goff[i]=lsize*i+headerSize;
	}

	goff[numproc]=fsize;

	for(i=1;i<numproc;i++)
	{
		current_line =(char*)calloc(1000,sizeof(char));
		MPI_File_read_at_all(mpi_filed, (MPI_Offset)goff[i], current_line, 1000, MPI_CHAR, &status);
		j=0;
		while(j<fsize && current_line[j] != '\n'){
			j++;
		}
		goff[i]+=(j+1);
		free(current_line);
	}

	return goff;
}

void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads){

	char *currentCarac;
	char * currentLine =(char*)calloc(MAX_LINE_SIZE,sizeof(char));
	int treshold2 = 0;
	int quality;
	unsigned int i, chr, nbchr = 0, mchr;
	int lastChr = -1;
	char forward;
	int next;
	size_t lineSize, offset_read_in_source_file;
	size_t size_data = 0;
	size_t size_read = 0;
	size_t coord, mcoord;
	size_t *readNumberByChr;
	size_t counter = 0;
	Read **reads = *preads;
	size_t bad_quality = 0;

	//we take the first line *
	//before calling parser paired, we know that localdata is at the begining of a read

	size_data = strlen(localData);

	next = tokenizer(localData,'\n', currentLine);

	assert(currentLine != NULL);

	offset_read_in_source_file = start_offset;


	nbchr = nbchrom;
	readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

	while(next && size_read < size_data){

		lineSize = strlen(currentLine) + 1;
		size_read += strlen(currentLine);

		assert(lineSize != 0);

		//we update the offset in the
		//source file
		currentLine[lineSize - 1] = '\n';
		currentLine[lineSize] = '\0';


		//Getting forward

		//GO TO FLAG
		currentCarac = strstr(currentLine, "\t");
		//GO TO RNAME (Chr name)
		currentCarac = strstr(currentCarac+1, "\t");
		if(lastChr == (nbchr - 1) && forward==0)
		{
			chr = (nbchr -1);
		}
		else
		{
			chr = getChr(currentCarac, chrNames, nbchr);
		}

		//GO TO COORD
		currentCarac = strstr(currentCarac+1, "\t");
		//TAKE COORD AND GO TO MAPQ
		coord = strtoull(currentCarac, &currentCarac, 10);

		//TAKE MAPQ AND GO TO CIGAR
		quality = strtoull(currentCarac, &currentCarac, 10);

		//GO TO RNEXT
		currentCarac = strstr(currentCarac+1, "\t");
		if(currentCarac[1]=='='){
			mchr = chr;
		}
		else if(currentCarac[1] == '*'){
			mchr = (nbchr-1);
		}
		else{
			mchr = getChr(currentCarac, chrNames, nbchr);
		}

		//GO TO PNEXT
		currentCarac = strstr(currentCarac+1, "\t");
		if(currentCarac == NULL){
			printf("nbread parsed : %zu / offset %zu\n",counter, offset_read_in_source_file);
			printf("mcoord : %s",currentLine);
			MPI_Abort(MPI_COMM_WORLD,77);
		}
		mcoord = strtoull(currentCarac, &currentCarac, 10);

		//forward values:
		//forward : 1
		//backward : 0
		//duplicate : 2
		forward = ((coord<mcoord)) ? 1 : ( (coord==mcoord) ? 2 : 0);

		if (chr < nbchr-1){

			if (quality > treshold2){
				reads[chr]->next = (Read*)malloc(sizeof(Read));
				assert(reads[chr]->next);
				//Flags
				reads[chr]->next->flags.is_mate=((forward==0)? 1 : 0);
				reads[chr]->next->flags.replace_gene_with_mgene=0;
				reads[chr]->next->flags.chr=chr;
				reads[chr]->next->flags.left=0;
				//Read
				reads[chr]->next->coord = coord;
				reads[chr]->next->quality = quality;
				//Mate
				reads[chr]->next->mcoord = mcoord;
				reads[chr]->next->mchr=mchr;
				reads[chr]->next->gene = 0;
				reads[chr]->next->mgene = 0;
				//File
				reads[chr]->next->offset_source_file=offset_read_in_source_file;
				reads[chr]->next->offset = lineSize;
				reads[chr] = reads[chr]->next;
				reads[chr]->next = NULL;
				readNumberByChr[chr]++;
			}

			else{
				bad_quality++;
			}
		}

		else{
			reads[nbchr-1]->next = (Read*)malloc(sizeof(Read));
			reads[nbchr-1]->next->flags.is_mate=0;
			reads[nbchr-1]->next->flags.replace_gene_with_mgene=0;
			reads[nbchr-1]->next->flags.chr=0;
			reads[nbchr-1]->next->flags.left=0;
			reads[nbchr-1]->next->coord = coord;
			reads[nbchr-1]->next->mcoord = mcoord;
			reads[nbchr-1]->next->gene = 0;
			reads[nbchr-1]->next->mgene = 0;
			reads[nbchr-1]->next->quality = 0;
			reads[nbchr-1]->next->mchr = 0;
			reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
			reads[nbchr-1]->next->offset = lineSize;
			reads[nbchr-1] = reads[nbchr-1]->next;
			reads[nbchr-1]->next = NULL;
			readNumberByChr[nbchr-1]++;
		}

		//we update the offset_read_in_source_file
		offset_read_in_source_file += lineSize;
		//we read the next line

		for(i=0;i<strlen(currentLine);i++){
			currentLine[i]=0;
		}
		next = tokenizer(NULL, '\n', currentLine);
		counter++;
	}

	for(i=0;i<nbchr;i++){
		preadNumberByChr[0][i] += readNumberByChr[i];
	}

	free(currentLine);
	free(readNumberByChr);

}

void parser_single(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads){

		char *currentCarac;
		char currentLine[MAX_LINE_SIZE];
		unsigned char quality;
		unsigned int i, chr, nbchr = 0;
		int lastChr = -1;
		int next;
		size_t lineSize, offset_read_in_source_file;
		size_t coord;
		size_t *readNumberByChr;
		size_t counter = 0;
		Read **reads = *preads;


		for(i=0;i<MAX_LINE_SIZE;i++){
			currentLine[i]=0;
		}

		//we take the first line *
		//before calling parsepaired, we know that localdata is at the begining of a read

		next = tokenizer(localData,'\n', currentLine);

		offset_read_in_source_file = start_offset;

		nbchr = nbchrom;
		readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

		while(next){

			lineSize = strlen(currentLine) + 1;

			//we update the offset in the
			//source file

			currentLine[lineSize - 1] = '\n';
			currentLine[lineSize] = '\0';


			//GO TO FLAG
			currentCarac = strstr(currentLine, "\t");

			//GO TO RNAME (Chr name)
			currentCarac = strstr(currentCarac+1, "\t");
			if(lastChr == (nbchr - 1))
			{
				chr = (nbchr -1);
			}
			else
			{
				chr = getChr(currentCarac, chrNames, nbchr);
			}


			//GO TO COORD
			currentCarac = strstr(currentCarac+1, "\t");
			//TAKE COORD AND GO TO MAPQ
			coord = strtoull(currentCarac, &currentCarac, 10);

			//TAKE MAPQ AND GO TO CIGAR
			quality = strtoull(currentCarac, &currentCarac, 10);

			if (chr < nbchr-1){

				if(quality > threshold){
					reads[chr]->next = NULL;
					reads[chr]->next = (Read*)malloc(sizeof(Read));
					assert(reads[chr]->next);
					//Flags
					reads[chr]->next->flags.is_mate=0;
					reads[chr]->next->flags.replace_gene_with_mgene=0;
					reads[chr]->next->flags.chr=chr;
					reads[chr]->next->flags.left=0;
					//Read
					reads[chr]->next->coord = coord;
					reads[chr]->next->quality = quality;
					//Mate
					reads[chr]->next->mcoord = 0;
					reads[chr]->next->mchr=0;
					reads[chr]->next->gene = 0;
					reads[chr]->next->mgene = 0;
					//File
					reads[chr]->next->offset_source_file=offset_read_in_source_file;
					reads[chr]->next->offset = lineSize;
					reads[chr]->next->next=NULL;
					reads[chr] = reads[chr]->next;
					readNumberByChr[chr]++;

				}
			}

			else{
				reads[nbchr-1]->next = (Read*)malloc(sizeof(Read));
				reads[nbchr-1]->next->flags.is_mate=0;
				reads[nbchr-1]->next->flags.replace_gene_with_mgene=0;
				reads[nbchr-1]->next->flags.chr=0;
				reads[nbchr-1]->next->flags.left=0;
				reads[nbchr-1]->next->coord = 0;
				reads[nbchr-1]->next->mcoord = 0;
				reads[nbchr-1]->next->gene = 0;
				reads[nbchr-1]->next->mgene = 0;
				reads[nbchr-1]->next->quality = 0;
				reads[nbchr-1]->next->mchr = 0;
				reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
				reads[nbchr-1]->next->offset = lineSize;
				reads[nbchr-1]->next->next = NULL;
				reads[nbchr-1] = reads[nbchr-1]->next;
				readNumberByChr[nbchr-1]++;
			}

			//we update the offset_read_in_source_file
			offset_read_in_source_file += lineSize;
			//we read the next line

			for(i=0;i<MAX_LINE_SIZE;i++){
				currentLine[i]=0;
			}
			next = tokenizer(NULL, '\n', currentLine);

			counter++;
		}


		for(i=0;i<nbchr;i++){
			preadNumberByChr[0][i] += readNumberByChr[i];
		}


}

int getChr(char *str, char** chrNames, int nbchr){
	int i=0, found=0, size;
	char *str1 = str, *str2;

	str2 = str1+1;

	for(; *str2 != '\t'; str2++);

	size=strlen(str1)-strlen(str2);

	char * tmp_chr =(char*)malloc(sizeof(char)*size);

	for(i=0;i<size;i++){
		tmp_chr[i]=str1[i+1];
	}
	tmp_chr[size-1]=0;

	assert(strlen(tmp_chr) != 0);
	for(i = 0, found = 0; i < nbchr && !found; i++){
		found = !strcmp(tmp_chr, chrNames[i]);
	}

	free(tmp_chr);
	return i-1;
}
