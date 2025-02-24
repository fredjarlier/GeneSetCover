/**
 * \file mergeSort.c
 * \author Nicolas Fedy to line 57. Then Thomas MAGALHAES, Paul PANGANIBAN, Nizar HDADECH
 * \date August 26th 2015
 */
#include "mergeSort.h"


Read* mergeSort(Read* c, size_t n){
	size_t q, p;
	Read* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSort(c, p);
		if(q > 1)
			mergeSort(d, q);
	}
	else
		d = c->next;
	d = structMerge(c, p, d, q);

	return d;
}

Read* structMerge(Read* c, size_t p, Read* d, size_t q){

	Read* t;

	while(1){
		if (c->next->coord > d->next->coord || (c->next->coord == d->next->coord && c->next->mcoord > d->next->mcoord) ){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}


Read* mergeSort_mcoord(Read* c, size_t n){
	size_t q, p;
	Read* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSort_mcoord(c, p);
		if(q > 1)
			mergeSort_mcoord(d, q);
	}
	else
		d = c->next;
	d = structMerge_mcoord(c, p, d, q);

	return d;
}

Read* structMerge_mcoord(Read* c, size_t p, Read* d, size_t q){

	Read* t;

	while(1){
		if (c->next->mcoord > d->next->mcoord || (c->next->mcoord == d->next->mcoord && c->next->coord > d->next->coord) ){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}

Read_chain* mergeSortReadChain(Read_chain* c, size_t n){
	size_t q, p;
	Read_chain* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortReadChain(c, p);
		if(q > 1)
			mergeSortReadChain(d, q);
	}
	else
		d = c->next;
	d = structMergeReadChain(c, p, d, q);

	return d;
}

Read_chain* structMergeReadChain(Read_chain* c, size_t p, Read_chain* d, size_t q){
	Read_chain* t;

	while(c->next && d->next && d->next->next){
		if (c->next->reads->mcoord > d->next->reads->mcoord){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}

Read* mergeSortGene(Read* c, size_t n){
	size_t q, p;
	Read* d;
	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortGene(c, p);
		if(q > 1)
			mergeSortGene(d, q);
	}
	else
		d = c->next;
	d = structMergeGene(c, p, d, q);

	return d;
}

Read* mergeSortMgene(Read* c, size_t n){
	size_t q, p;
	Read* d;
	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortMgene(c, p);
		if(q > 1)
			mergeSortMgene(d, q);
	}
	else
		d = c->next;

	d = structMergeMgene(c, p, d, q);

	return d;
}
Read* structMergeGene(Read* c, size_t p, Read* d, size_t q){

	Read* t;

	while(c->next && d->next){


		if ( c->next->gene > d->next->gene|| (c->next->gene == d->next->gene && c->next->mgene > d->next->mgene) ){

			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}
Read* structMergeMgene(Read* c, size_t p, Read* d, size_t q){

	Read* t;

	while(c->next && d->next){
		if (c->next->gene==d->next->gene && c->next->mgene > d->next->mgene){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}


Cluster* mergeSortCluster_by_pnbread_inverted(Cluster* c, size_t n){
	size_t q, p;
	Cluster* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortCluster_by_pnbread_inverted(c, p);
		if(q > 1)
			mergeSortCluster_by_pnbread_inverted(d, q);
	}
	else
		d = c->next;
	d = structMergeCluster_by_pnbread_inverted(c, p, d, q);

	return d;
}

Cluster* structMergeCluster_by_pnbread_inverted(Cluster* c, size_t p, Cluster* d, size_t q){

	Cluster* t;

	while(1){


		if (c->pnbreads < d->pnbreads){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}


Cluster* mergeSortCluster_by_Mgene(Cluster* c, size_t n){
	size_t q, p;
	Cluster* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortCluster_by_Mgene(c, p);
		if(q > 1)
			mergeSortCluster_by_Mgene(d, q);
	}
	else
		d = c->next;
	d = structMergeCluster_by_Mgene(c, p, d, q);

	return d;
}

Cluster* structMergeCluster_by_Mgene(Cluster* c, size_t p, Cluster* d, size_t q){

	Cluster* t;

	while(1){
		if (c->next->read->mgene > d->next->read->mgene ){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}



Cluster* mergeSortCluster_by_Gene(Cluster* c, size_t n){
	size_t q, p;
	Cluster* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortCluster_by_Gene(c, p);
		if(q > 1)
			mergeSortCluster_by_Gene(d, q);
	}
	else
		d = c->next;
	d = structMergeCluster_by_Gene(c, p, d, q);

	return d;
}

Cluster* structMergeCluster_by_Gene(Cluster* c, size_t p, Cluster* d, size_t q){

	Cluster* t;

	while(1){
		if ((c->next->read->gene > d->next->read->gene)
		|| (c->next->read->gene == d->next->read->gene && c->next->read->mgene > d->next->read->mgene )){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}




ClusterSplit* mergeSortClusterSplit(ClusterSplit* c, size_t n){
	size_t q, p;
	ClusterSplit* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortClusterSplit(c, p);
		if(q > 1)
			mergeSortClusterSplit(d, q);
	}
	else
		d = c->next;
	d = structMergeClusterSplit(c, p, d, q);

	return d;
}

ClusterSplit* structMergeClusterSplit(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q){

	ClusterSplit* t;

	while(1){


		if (c->next->nbread > d->next->nbread){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}



ClusterSplit* mergeSortClusterSplit_inverted(ClusterSplit* c, size_t n){
	size_t q, p;
	ClusterSplit* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortClusterSplit_inverted(c, p);
		if(q > 1)
			mergeSortClusterSplit_inverted(d, q);
	}
	else
		d = c->next;
	d = structMergeClusterSplit_inverted(c, p, d, q);

	return d;
}

ClusterSplit* structMergeClusterSplit_inverted(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q){

	ClusterSplit* t;

	while(1){


		if (c->next->nbread < d->next->nbread){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}

ClusterSplit* mergeSortClusterSplit_by_gene(ClusterSplit* c, size_t n){
	size_t q, p;
	ClusterSplit* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortClusterSplit_by_gene(c, p);
		if(q > 1)
			mergeSortClusterSplit_by_gene(d, q);
	}
	else
		d = c->next;
	d = structMergeClusterSplit_by_gene(c, p, d, q);

	return d;
}

ClusterSplit* structMergeClusterSplit_by_gene(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q){

	ClusterSplit* t;

	while(1){


		if (c->next->gene > d->next->gene
			|| (c->next->gene == d->next->gene  && c->next->mgene > d->next->mgene)){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}
ClusterSplit* mergeSortClusterSplit_by_rank(ClusterSplit* c, size_t n){
	size_t q, p;
	ClusterSplit* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSortClusterSplit_by_rank(c, p);
		if(q > 1)
			mergeSortClusterSplit_by_rank(d, q);
	}
	else
		d = c->next;
	d = structMergeClusterSplit_by_rank(c, p, d, q);

	return d;
}

ClusterSplit* structMergeClusterSplit_by_rank(ClusterSplit* c, size_t p, ClusterSplit* d, size_t q){

	ClusterSplit* t;

	while(1){


		if (c->next->maxrank > d->next->maxrank){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}//if
			--p;
		}//else
		c = c->next;
	}//while(1)
	return d;
}
Iteration* mergeSort_iteration(Iteration* c, size_t n){
	size_t q, p;
	Iteration* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSort_iteration(c, p);
		if(q > 1)
			mergeSort_iteration(d, q);
	}
	else
		d = c->next;
	d = structMerge_iteration(c, p, d, q);

	return d;
}

Iteration* structMerge_iteration(Iteration* c, size_t p, Iteration* d, size_t q){

	Iteration* t;

	while(1){
		if (c->rate > d->rate){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}

void merge_coord (size_t *a, size_t n, size_t m) {
	size_t i, j, k;
    size_t *x = malloc(n * sizeof (size_t));
    for (i = 0, j = m, k = 0; k < n; k++) {
        x[k] = j == n      ? a[i++]
             : i == m      ? a[j++]
             : a[j] < a[i] ? a[j++]
             :               a[i++];
    }
    for (i = 0; i < n; i++) {
        a[i] = x[i];
    }
    free(x);
}

void merge_sort_coord (size_t *a, size_t n) {
    if (n < 2)
        return;
    size_t m = n / 2;
    merge_sort_coord(a, m);
    merge_sort_coord(a + m, n - m);
    merge_coord(a, n, m);
}
