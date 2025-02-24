/**
 * \file duplicata.c
 * \author Paul PANGANIBAN
 * \date August 26th 2015
 */

#include "duplicata.h"


size_t duplicata(Cluster * clust){

	Cluster * cluster = clust;
	Read * c;
	Read* d = NULL;

	size_t diff = 0;
	size_t result = 0;
	while(cluster){

		c = cluster->read;
		diff=0;
		while(c->next){

				d = c->next ;
				//if DUPLICATA, we delete next (read d)
				if(c->coord == d->coord
					&& c->mcoord != 0
					&& d->mcoord != 0
					&& c->mcoord == d->mcoord
					){
					c->next = d->next;
					free(d);
					cluster -> pnbreads --;
					diff++;

				}//if duplicata

				//else we look next
				else{
					c = c->next;
				}
		}
		cluster = cluster ->next;
		result +=diff;
	}

	return result;

}
