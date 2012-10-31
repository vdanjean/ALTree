#ifndef _STATS_H
#define _STATS_H

#include "datatype.h"

#define SAMPLESIZE 5
#define Seuil_ONLY_CASE 3
#define PERM 1000

struct cc {
	datatype_t cases;
	datatype_t controls;
};

struct classical_chi2_res {
	datatype_t chi2;
	int chi2invalid;
	int error;
	int sum_control;
	int sum_case;
};

struct classical_chi2_res classical_chi2(int nb_nodes, struct cc *nodes);
void definition_p_chi2(datatype_t p, datatype_t pprop);
int chi2_significatif(int ddl, datatype_t chi2);
int chi2_fisher_significatif(datatype_t pvalue);
datatype_t reech_chi2(int sum_case, int sum_control,
		      int nb_nodes, int chi2_reel, struct cc *nodes);
int reech_significatif(datatype_t p_val);
void random_clades(int nb_nodes, struct cc *nodes,
		   int cases, int controls, struct cc *clades);

#endif

