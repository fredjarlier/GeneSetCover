/**
* \file tree.c
* \date 27 aout 2015
* \author Paul PANGANIBAN
*/
#include "tree.h"

int father(int rank){
	
	if( rank%2 == 0 ){
		return (rank/2 -1);
	}
	else{
		return (rank/2); 
	}
}

void child(int rank, size2_t * child, int nbProcess){
	child->a =0;
	child->b =0;
	if (rank*2+1<nbProcess)
		child -> a = (size_t)rank*2+1;
	if (rank*2+2<nbProcess)
		child -> b = (size_t)rank*2+2;
}

int level(int rank){
	int level = -1 ;
	double i = 0;
	
	if (rank == 1){
		return 1;	
	}
	else {
		while (level == -1){
			double p = pow(2,i);
			if(rank <= p ){
				level = i;
			}

			else{
				i++;
			}
		}
		return level; 
	}

}
