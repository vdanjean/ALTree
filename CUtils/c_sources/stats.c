
#include "debug.h"
#include "stats.h"
#include "fisher.h"
#include "chisq.h"
#include "myrand.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>
#include <strings.h>

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

static datatype_t *chi2_seuil=NULL;
static int nb_seuils=0;

static datatype_t chi2_p=-1;
static datatype_t test_prop_p=-1;

//Utilisé aussi pour le test F
void definition_p_chi2(datatype_t p, datatype_t pprop) {
	if (p != -1) {
		chi2_p=p;
	}
	if (pprop != -1) {
		test_prop_p=pprop;
	}
}

int chi2_significatif(int ddl, datatype_t chi2) {
	if (ddl < 1) {
		fprintf(stderr, "Warning: chi[%i] asked...\n", ddl);
	}
	if (ddl >= nb_seuils) {
		chi2_seuil = realloc(chi2_seuil, (ddl+1)*sizeof(datatype_t));
		bzero(chi2_seuil+nb_seuils, (ddl+1-nb_seuils)*sizeof(datatype_t));
		nb_seuils=ddl+1;
	}
	if (chi2_seuil[ddl]==0.0) {
		assert(chi2_p >= 0.0);
		chi2_seuil[ddl]=critchi(chi2_p, ddl);
	}
	return (chi2 > chi2_seuil[ddl]);
}

int chi2_fisher_significatif(datatype_t pvalue)
// "Meme fonction pour le test F, je ne la ré-écris pas..." !?!?
{
	return pvalue < chi2_p;
}

datatype_t reech_chi2(int sum_case, int sum_control,
		      int nb_nodes, int chi2_reel, struct cc *nodes)
{
	int sum_total=sum_case+sum_control;
	datatype_t p_val=0.0;

	datatype_t th_c[nb_nodes];
	datatype_t th_m[nb_nodes];

	int i,k;
	for(i=0; i<nb_nodes; i++) {
		th_c[i]=((datatype_t)(sum_control*(nodes[i].cases+nodes[i].controls)))/sum_total;
		th_m[i]=((datatype_t)(sum_case*(nodes[i].cases+nodes[i].controls)))/sum_total;
	}

	struct cc clades[nb_nodes];
	for (k=1;k<=PERM; k++){
		random_clades(nb_nodes, nodes, sum_case, sum_control, clades);

		datatype_t chi2=0.0;
		for(i=0; i<nb_nodes; i++) {
			datatype_t t;
			t=clades[i].cases - th_m[i];
			chi2 += t*t/th_m[i];
			t=clades[i].controls - th_c[i];
			chi2 += t*t/th_c[i];
		}
		if (chi2 >= chi2_reel){
			p_val++;
		}
	}
	//debug("pval=%g\n", p_val/PERM);
	return (p_val/PERM);
}

int reech_significatif(datatype_t p_val)
{
	return (p_val<=chi2_p);
}

void random_clades(int nb_nodes, struct cc *nodes,
		   int cases, int controls, struct cc *clades)
{
	bzero(clades, nb_nodes*sizeof(struct cc));
	int c;
	for(c=0; c<nb_nodes; c++) {
		int i;
		int s=nodes[c].cases+nodes[c].controls;
		for(i=0; i<s; i++) {
			int alea=myrand(cases+controls);
			if (alea < cases) {
				cases--;
				clades[c].cases++;
			} else {
				controls--;
				clades[c].controls++;
				}
		}
	}
}
