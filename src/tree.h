/**
 * \file tree.h
 * \date 27 aout 2015
 * \brief Regroups all methods related to creating a binary tree
 * \author Paul PANGANIBAN
 */
#ifndef TREE_H
#define TREE_H

#include <math.h>
#include "genes.h"

/**
 * \brief Tells a proc which rank is his father
 *
 * \param rank rank of process
 * \return rank of the father
 */
int father(int rank);

/**
 * \brief Tells a proc which ranks are his children
 *
 * \param[in] rank rank of process
 * \param[out] child structure which will save the two children
 * \param[in] nbProcess number of process
 */
void child(int rank, size2_t * child, int nbProce);

/**
 * \brief Tells a proc which level he has in the binary tree.
 * 		0 is the tree root.
 *
 * \param rank rank of process
 * \return level of the process
 */
int level(int rank);

#endif
