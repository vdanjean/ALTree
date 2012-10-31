
#include "debug.h"
#include "stats.h"
#include <stdio.h>
#include <stdlib.h>

struct classical_chi2_res classical_chi2(int nb_nodes, struct cc *nodes) {
	struct classical_chi2_res res={
		.chi2=0.0,
		.chi2invalid=0,
		.error=0,
		.sum_control=0,
		.sum_case=0,
	};

	int i;
	for (i=0; i<nb_nodes; i++) {
		res.sum_control+=nodes[i].controls;
		res.sum_case+=nodes[i].cases;
	}
	int sum_total=res.sum_control+res.sum_case;
	int ddl=nb_nodes-1;

	if (ddl==0) { // 1 seul clade
		res.error=4;
		return res;
	}
	if (res.sum_case==0) {
		res.error=1;
		return res;
	}
	if (res.sum_control==0) {
		res.error=2;
		return res;
	}
	
	for(i=0; i<nb_nodes; i++) {
		int m=nodes[i].cases;
		int c=nodes[i].controls;
		if (m==0 && c==0) {
			fprintf(stderr, "no case and no control in a node!\n");
			exit(1);
		}
		datatype_t t_m=((double)((m+c)*res.sum_case))/sum_total;
		res.chi2 += (m - t_m)*(m - t_m)/t_m;

		datatype_t t_c=((double)((m+c)*res.sum_control))/sum_total;
		res.chi2 += (c - t_c)*(c - t_c)/t_c;

		if( t_m <= SAMPLESIZE
		    || t_c <= SAMPLESIZE ) {
			res.chi2invalid++;
		}
	}
	return res;
}

