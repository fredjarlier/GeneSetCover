/**
 * \file parser.h
 * \brief Contains all functions for converting line file in struct Read
 * \author Nicolas FEDY, Leonor SIROTTI, Paul PANGANIBAN, Thomas MAGALHAES, Nizar HDADECH
 * \date August 26th 2015
 */

#ifndef PARSER
	#define PARSER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"
#include "tokenizer.h"
#include "time.h"

#define MAX_LINE_SIZE 2048
#define UNMAPPED "unmapped"

typedef struct Flags Flags;
typedef struct Read Read;

struct Flags
{
	unsigned char chr : 5;
	unsigned char is_mate : 1;
	unsigned char left : 1;
	unsigned char replace_gene_with_mgene:1;
};

struct Read
{
	size_t offset_source_file;
	size_t offset;

	size_t coord;
	size_t gene;

	size_t mcoord;
	size_t mgene;

	Flags flags;

	char quality;
	char mchr;

	char exon_coord;
	char exon_mcoord;
	Read* next;
	Read* link;
}; //attribute packed slows by x2 times


typedef struct Read_chain
{
	struct Read* reads;
	struct Read_chain* next;
}Read_chain;

/**
 * \brief Find the header, its size and the chromosomes name
 *
 * \fn unsigned int find_header(char *localData, int rank, size_t *unmappedSize, int *pnbchr, char **pheader, char ***pchrNames)
 * \param localData The local data read directly from the input file
 * \param rank Rank of the process
 * \param[out] unmappedSize Size of lines that are not "Read"
 * \param[out] pnbchr The number of chromosome. Will be set in this function.
 * \param[out] pheader Header in char *
 * \param[out] pchrNames Reference array to the chromosomes names. Will be filled during the function.
 *
 * \return Size of the header, chromosome name
 */
unsigned int find_header(char *localData, int rank, size_t *unmappedSize, int *pnbchr, char **pheader, char ***pchrNames);

/**
 * \brief Initialize the start offset for each process
 *
 * \fn size_t * init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank)
 * \param mpi_filed File to read
 * \param headerSize The size of the header
 * \param fsize The size of the mpi_filed
 * \param numproc The number of process
 * \param rank Rank of the process
 *
 * \return An array of "numproc" cells with the first offset to read for each process (last cell is the last offset of mpi_filed = fsize)
 */
size_t * init_goff(MPI_File mpi_filed,unsigned int headerSize,size_t fsize,int numproc,int rank);

/**
 * \brief Parse localData to reads. This functions takes care of the mate of each read.
 *
 * \fn void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads)
 * \param localData The local data read directly from the input file
 * \param rank Rank of the process
 * \param start_offset First localData offset in the file
 * \param threshold Minimum quality that the reads have to be. If the read quality isn't high enough, it won't be keeped.
 * \param nbchrom The number of chromosome
 * \param[out] preadNumberByChr Reference array to the number of reads by chromosome.
 * \param chrNames Reference array to the chromosomes names.
 * \param[out] preads The reads linked list. Will be set in this function.
 *
 * \return void
 */
void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads);

/**
 * \brief Parse localData to reads. This functions doesn't take care of the mate of each read.
 *
 * \fn void parser_single(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads)
 * \param localData The local data read directly from the input file
 * \param rank Rank of the process
 * \param start_offset First localData offset in the file
 * \param threshold Minimum quality that the reads have to be. If the read quality isn't high enough, it won't be keeped.
 * \param nbchrom The number of chromosome
 * \param[out] preadNumberByChr Reference array to the number of reads by chromosome.
 * \param chrNames Reference array to the chromosomes names.
 * \param[out] preads The reads linked list. Will be set in this function.
 *
 * \return void
 */
void parser_single(char *localData, int rank, size_t start_offset, unsigned char threshold,int nbchrom, size_t **preadNumberByChr, char ** chrNames, Read ***preads);

/**
 * \brief Extract integer number of current chromosome given as a string
 *
 * \param str The string containing the chromosome ID
 * \param chrNames Reference array to the chromosomes names
 * \param nbchr The number of chromosome
 *
 * \return The chromosome as an integer
 */
int getChr(char* str, char** chrNames, int nbchr);

#endif
