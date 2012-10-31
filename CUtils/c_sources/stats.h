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

struct calcul_chi2_res {
	datatype_t chi2;
	datatype_t p_val;
	int error;
	int significatif;
	char *texte;
	char *warning;
};

struct classical_chi2_res classical_chi2(int nb_nodes, const struct cc *nodes);
struct calcul_chi2_res calcul_chi2(int nb_nodes, const struct cc *nodes,
				   int sign_util, int texte, struct cc *th);
void definition_p_chi2(datatype_t p, datatype_t pprop);
int chi2_significatif(int ddl, datatype_t chi2);
int chi2_fisher_significatif(datatype_t pvalue);
/* nodes and th: array of size nb_nodes.
   nodes is input only, th is used internally
*/
datatype_t reech_chi2(int sum_case, int sum_control, int nb_nodes,
		      int chi2_reel, const struct cc *nodes, struct cc *th);
int reech_significatif(datatype_t p_val);
void random_clades(int nb_nodes, const struct cc *nodes,
		   int cases, int controls, struct cc *clades);

#endif

