/**
 * \file genes.h
 * \brief All functions to find gene of a read and treatment on gene file
 * \author Thomas MAGALHAES, Paul PANGANIBAN
 * \date August 26th 2015
 */

#ifndef GENES
#define GENES

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "parser.h"

typedef struct size2_t
{
	size_t a;
	size_t b;
}size2_t;



/**
 * \brief Set gene for each read in "Read* read" with buffer_fil.
 * 		  Delete reads that didnt match
 * \param buffer_file Char * with gene ID
 * \param[out] backward Linked list of Read that we will treat
 * \return a struct*size2_t with attributes "a" equals to the number of reads which matched and attributes "b" the number of unmatched
 */
size2_t getGeneAsNumber(char* buffer_file, Read* read);


/**
 * \brief Set mate gene for each backward in "Read_chain * backward" with buffer_file
 * \param buffer_file Char * with gene ID
 * \param[out] backward Linked list of Read_chain that we will treat
 * \return a struct*size2_t with attributes "a" equals to the number of backwards which matched and attributes "b" the number of unmatched
 */
size2_t getMGeneAsNumber(char* buffer_file, Read_chain* backward);

/**
 * \brief Find bounds of gene in char * line
 * \param line Contains the string from which the gene number should be extracted
 * \param[out] size Will contain the start and end bounds of the gene found.
 * \return the gene value as a number
 */
size_t getGeneLineBounds(char* line, size2_t* size);

/**
 * \brief Gets every backward read from @reads
 * \param[out] reads The read array containing the entire reads
 * \param nbchr The number of chromosomes "reads" is made of
 * \param[out] count An array for counting the backwards / chromosomes
 * \return A linked list of Read_chain of each backwards by chromosomes.
 */
Read_chain** getBackward(Read** reads, int nbchr, size_t* count);

/**
 * \brief Free each chromosome's backward list
 * \param[out] chain the array of reads indexed by chr
 * \param nbchr the number of chromosomes
 * \return void
 */
void freeBackward(Read_chain** chain, int nbchr);

#endif
