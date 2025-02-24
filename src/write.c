/**
 * \file write.c
 * \date 27 aout 2015
 * \author	Thomas Magalhaes
 * 			Paul Panganiban
 */

#include "write.h"

void getFusionsBoundaries(size_t* forward, size_t* backward, Cluster* c, size_t count_cluster)
{
	Cluster* tmp = c;
	size_t m;

	//Loop through clusters
	for(m = 0; m<count_cluster; m++)
	{
		forward[m]=0; //Set it to the lower
		backward[m]=0-1; //Set it to the bigger
		Read* r = tmp->read;
		while(r)
		{
			if(!r->flags.is_mate)
			{
				if(r->coord > forward[m])
					forward[m] = r->coord;
				if(r->mcoord < backward[m])
					backward[m] = r->mcoord;
			}
			else
			{
				if(r->coord < backward[m])
					backward[m] = r->coord;
				if(r->mcoord > forward[m])
					forward[m] = r->mcoord;
			}

			r = r->next;
		}

		tmp = tmp->next;
	}
}

void writeCluster(int rank, int num_proc, Cluster* clusters, size_t count_clusters, char** chrNames, char* fileName, char* header, char* outputDir, char* exon_valid_cluster)
{
	int i;
	MPI_Status status;
	Cluster* tmp;
	size_t recv_offs[num_proc];

	//Bcast header length for each process to know where does the reads start
	int headerLength = 0;
	if(!rank)
		headerLength=strlen(header);
	MPI_Bcast(&headerLength, 1, MPI_INT, 0, MPI_COMM_WORLD);

	char *data = NULL;
	//Read reads
	size_t *clSize = calloc(count_clusters, sizeof(size_t));
	clSize[0] = 0;

	//Read the data
	size_t dataSize = writeCluster_read(rank, num_proc, clusters, &data, clSize, count_clusters, fileName);

	tmp = clusters;

	//Define for each process where it has to start writing its reads
	size_t offset_count_reads=0;
	MPI_Gather(&dataSize, 1, MPI_UNSIGNED_LONG, recv_offs, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if(!rank)
	{
		size_t rt[num_proc];
		rt[0] = headerLength;
		for(i = 1; i<num_proc; i++)
		{
			rt[i] = rt[i-1] + recv_offs[i-1];
		}
		rt[0] = headerLength;
		for(i = 0; i<num_proc; i++)
		{
			recv_offs[i] = rt[i];
		}
	}
	MPI_Scatter(recv_offs, 1, MPI_UNSIGNED_LONG, &offset_count_reads, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

	//Cluster
	size_t forwards[count_clusters];
	size_t backwards[count_clusters];
	getFusionsBoundaries(forwards, backwards, clusters, count_clusters);
	char** cW = malloc(sizeof(char*)*count_clusters);
	size_t offset_count = 0;
	size_t count_tmp = offset_count_reads;

	//Load clusters strings
	tmp = clusters;
	for(i = 0; i < count_clusters; i++)
	{
		cW[i] = malloc(255);
		//ChrA	ChrB	GeneA(ID)	GeneB(ID)	nbreads	offset_start	offset_length	FBStart	FBEnd	Exon_match_percent
		sprintf(cW[i], "%s\t%s\t%s(%zu)\t%s(%zu)\t%zu\t%zu\t%zu\t%zu\t%zu\t%d\n", chrNames[(int)tmp->read->flags.chr], chrNames[(int)tmp->read->mchr], tmp->gene_name,tmp->read->gene, tmp->mgene_name,tmp->read->mgene, tmp->pnbreads, count_tmp, clSize[i], forwards[i], backwards[i], (int)exon_valid_cluster[i]);
		count_tmp+=clSize[i];
		offset_count+=strlen(cW[i]);
		tmp = tmp->next;
	}
	free(clSize);

	//Define for each process where it has to write its clusters
	MPI_Gather(&offset_count, 1, MPI_UNSIGNED_LONG, recv_offs, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if(!rank)
	{
		size_t rt[num_proc];
		rt[0] = 0;
		for(i = 1; i<num_proc; i++)
		{
			rt[i] = rt[i-1] + recv_offs[i-1];
		}
		for(i = 0; i<num_proc; i++)
		{
			recv_offs[i] = rt[i];
		}
	}
	MPI_Scatter(recv_offs, 1, MPI_UNSIGNED_LONG, &offset_count, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

	//Write clusters
	MPI_File clusterFile;
	char clusterFileName[255];
	sprintf(clusterFileName, "%s/cluster_list.cl", outputDir);
	MPI_File_open(MPI_COMM_WORLD, clusterFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &clusterFile);
	MPI_Offset offset = (MPI_Offset)offset_count;
	MPI_File_seek(clusterFile, offset, MPI_SEEK_SET);
	i = 0;
	tmp = clusters;
	while(tmp)
	{
		MPI_File_write(clusterFile, cW[i], strlen(cW[i]), MPI_CHAR, &status);
		tmp = tmp->next;
		i++;
	}
	MPI_File_close(&clusterFile);

	//Free clusters data
	for(i = 0; i < count_clusters; i++)
	{
		free(cW[i]);
	}
	free(cW);

	//Write reads
	sprintf(clusterFileName, "%s/reads.sam", outputDir);
	MPI_File_open(MPI_COMM_WORLD, clusterFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &clusterFile);
	offset = (MPI_Offset)offset_count_reads;
	if(!rank)
	{
		MPI_File_seek(clusterFile, 0, MPI_SEEK_SET);
		MPI_File_write(clusterFile, header, strlen(header), MPI_CHAR, &status);
	}
	else
	{
		MPI_File_seek(clusterFile, offset, MPI_SEEK_SET);
	}
	MPI_File_write(clusterFile, data, dataSize, MPI_CHAR, &status);
	MPI_File_close(&clusterFile);

	//Free reads strings
	free(data);

	MPI_Barrier(MPI_COMM_WORLD);
}

size_t writeCluster_read(int rank, int num_proc, Cluster* clusters, char **data, size_t *clSize, size_t count_clusters, char* fileName)
{
	size_t count_size = 0;
	Cluster* tmp = NULL;
	Read *reads = NULL;
	int i = 0;
	int j = 0;

	//Define each cluster's size
	tmp = clusters;
	while(tmp != NULL)
	{
		reads = tmp->read;
		while(reads)
		{
			count_size += reads->offset;
			clSize[i] += reads->offset;
			reads = reads->next;
		}
		i++;
		tmp = tmp->next;
	}

	(*data) = (char *)calloc( (count_size+1), sizeof(char));

	//Read the reads
	MPI_File file;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	if(rank == 0)
		printf("Opened file %s\n", fileName);
	tmp = clusters;
	size_t count = 0;
	while(tmp)
	{
		reads = tmp->read;
		for(i = 0; i < tmp->pnbreads; i++)
		{
			MPI_File_read_at( file, (MPI_Offset)reads->offset_source_file,*(data)+count, reads->offset, MPI_CHAR, &status);
			count += reads->offset;
			(*data)[count]='\n';
			reads = reads->next;
		}
		tmp = tmp->next;
		j++;
	}
	MPI_File_close(&file);

	//Safety for the string to have an end
	data[0][count_size] = '\0';

	return count_size;
}
