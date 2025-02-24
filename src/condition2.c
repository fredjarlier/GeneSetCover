/**
 * \file condition2.c
 * \author Paul PANGANIBAN
 * \date August 26th 2015
 */
#include "condition2.h"



size_t dx(Read * c , Read *d){
	if(c->coord > d->coord){
		return (c->coord - d->coord);
	}

	else{
		return (d->coord - c->coord);
	}
}

size_t dy(Read * c , Read * d){
	if(c->mcoord > d->mcoord){
		return (c->mcoord - d->mcoord);
	}
	else{
		return (d->mcoord - c->mcoord);
	}

}

void condition2_on_clusters(Cluster ** pClusters , size_t alpha, size_t minimum_number_reads){


	Cluster * cl_curr = *pClusters;
	Cluster * cl_before = NULL;
	Read * forwards = NULL;
	Read * forwards_first = NULL;
	Read * backwards = NULL ;
	Read * backwards_first = NULL;
	Read * r_curr = NULL;

	size_t count_cl_supp = 0;
	size_t count_total_read_supp =0;
	size_t nb_forward = 0 ;
	size_t nb_backward = 0 ;

	Read * fake =(Read*)malloc(sizeof(Read));
	fake->coord = 0;
	fake->mcoord = 0;
	fake->next = NULL;
	//we will apply C2 on all clusters
	while(cl_curr){
		//we apply C2

		//we will seperate forward and backwards read
		forwards = NULL;
		forwards_first = NULL;
		backwards = NULL ;
		backwards_first = NULL;
		nb_backward = 0 ;
		nb_forward = 0;

		fake -> next = cl_curr->read;
		mergeSort(fake, cl_curr->pnbreads);
		cl_curr->read = fake ->next;
		r_curr=cl_curr->read;

		while(r_curr){
			if(r_curr->flags.is_mate == 0){
				if(!forwards){
					forwards = r_curr;
					forwards_first = r_curr;
				}
				else{
					forwards->next=r_curr;
					forwards = forwards->next;
				}
				nb_forward++;
			}
			else{
				if(!backwards){
					backwards = r_curr;
					backwards_first = r_curr;
				}
				else{
					backwards->next=r_curr;
					backwards =backwards->next;
				}
				nb_backward++;
			}
			r_curr=r_curr->next;
		}
		if(forwards){
			forwards->next = NULL;
		}
		if(backwards){
			backwards->next = NULL;
		}
		forwards = forwards_first;
		backwards = backwards_first;



		//C2 on forwards
		size_t nb_fw_supp = condtion2_on_read(&forwards,alpha);
		forwards_first=forwards;

		//C2 on backwards
		size_t nb_bw_supp = condtion2_on_read(&backwards,alpha);
		backwards_first=backwards;

		//we link forward and backwards
		//if there is no forwards cl_curr->read is only backwards
		if(!forwards){
			cl_curr->read=backwards_first;
		}
		else{
			//we find the last forwards
			while(forwards->next){
				forwards=forwards->next;
			}
			forwards->next=backwards_first;
			cl_curr->read=forwards_first;
		}

		count_total_read_supp+=nb_fw_supp;
		count_total_read_supp+=nb_bw_supp;

		cl_curr->pnbreads-=nb_fw_supp;
		cl_curr->pnbreads-=nb_bw_supp;


		//if this cluster have now less than 6 reads we delete it
		if(cl_curr->pnbreads <= minimum_number_reads){
			count_cl_supp++;
			//if its the first
			if(!cl_before){
				(*pClusters) = cl_curr->next;
				free(cl_curr);
				cl_curr=(*pClusters);
			}

			else{
				cl_before->next = cl_curr->next;
				free(cl_curr);
				cl_curr=cl_before->next;
			}
		}

		else{
			cl_before = cl_curr;
			cl_curr=cl_curr->next;
		}
	}//while cl_curr

	free(fake);
}

size_t condtion2_on_read(Read ** reads_cluster, size_t alpha){
	size_t nb_read_supp  = 0;
	Read * r_first = *reads_cluster;
	Read * r_curr = r_first;
	Read * r_before = NULL;
	Read * r_in_frame = NULL;

	size_t frame_pos_max = 0;
	size_t pos_rcurr = 0;
	size_t lmin = (1-alpha/2);
	size_t lmax = alpha/2;

	//we will check C2 for each reads
	while(r_curr){

		//we keep the pos of read curr
		pos_rcurr = r_curr->coord;
		frame_pos_max = pos_rcurr+lmax-lmin;

		//we place our read_frame
		if(r_curr->next){
			r_in_frame = r_curr->next;

			//we will check if r_curr verify C2 with all reads in the frame
			while(r_in_frame !=NULL && r_in_frame->coord <= frame_pos_max){

				//CONDITION 2  - we will use the bit left on attribute Flags in read to know if it pass C2
				//if they pass C2  we put the flag at 1
				if( (dx(r_curr,r_in_frame)+dy(r_curr,r_in_frame)) <= (lmax-lmin) ){
					r_curr->flags.left=1;
					r_in_frame->flags.left=1;
				}
				r_in_frame = r_in_frame->next;
			}
		}
		r_curr=r_curr->next;
	}

	//Now we will delete all read that doesnt pass C2
	r_curr = r_first;
	r_before = NULL;

	while(r_curr){
		//if it doesnt verify C2 we delete it
		if(r_curr->flags.left == 0){
			nb_read_supp++;
			//if its the first
			if(!r_before){
				*reads_cluster = r_curr->next;
				free(r_curr);
				r_curr = (*reads_cluster);
			}
			//if it's not the first
			else{
				r_before->next=r_curr->next;
				free(r_curr);
				r_curr=r_before->next;
			}
		}
		else{
			r_before = r_curr;
			r_curr = r_curr->next;
		}
	}

	return nb_read_supp;
}


