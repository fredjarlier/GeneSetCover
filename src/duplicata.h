/**
 * \file duplicata.h
 * \brief Contains duplicata function
 * \author  Paul PANGANIBAN
 * \date August 26th 2015
 */

#ifndef DUPLICATA_H
#define DUPLICATA_H

#include "stdio.h"
#include "stdlib.h"
#include "biSort.h"

/**
 * \brief Delete duplicates reads (2 reads with same coord and mcoord) in Clusters c
 * \param[out] c Linked list of Cluster
 * \return Number of reads deleted
 */
size_t duplicata(Cluster *c);

#endif
