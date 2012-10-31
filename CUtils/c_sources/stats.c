
#include "debug.h"
#include "stats.h"
#include "fisher.h"
#include "chisq.h"
#include "myrand.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>
#include <strings.h>
#include <string.h>

struct classical_chi2_res classical_chi2(int nb_nodes,
					 const struct cc *nodes) {
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

#define msprintf(str,...) \
	({						  \
		char* ch=NULL;						\
		size_t len=1+snprintf(ch, 0, str, ##__VA_ARGS__);	\
		ch=(char*)malloc(len);					\
		snprintf(ch, len, str, ##__VA_ARGS__);			\
		ch;							\
	})

#define msprintfcat(dest, str,...)					\
	({								\
		char* ch=NULL;						\
		size_t dlen=((dest)?strlen(dest):0);			\
		size_t slen=snprintf(ch, 0, str, ##__VA_ARGS__);	\
		ch=(char*)realloc(dest,dlen+slen+1);			\
		snprintf(ch+dlen, slen+1, str, ##__VA_ARGS__);		\
		ch;							\
	})

struct calcul_chi2_res calcul_chi2(int nb_nodes, const struct cc *nodes,
				   int sign_util, int texte, struct cc *th)
{
	struct classical_chi2_res r=
		classical_chi2(nb_nodes, nodes);

	struct calcul_chi2_res res={
		.error=r.error,
		.significatif=0,
		.texte=NULL,
		.warning=NULL,
	};

	assert(!(sign_util && !texte));
	if (res.error != 0) {
		if (!texte) { return res; }

		// TODO: A vérifier : est-ce OK de mettre $significatif à 0
		// la valeur est utilisée au retour de cette fonction
//		res.significatif=0;
		switch (res.error) {
		case 1:
			res.texte=msprintf("No cases,  (%i controls)", r.sum_control);
			break;
		case 2:
			res.texte=msprintf("No controls: only %i cases", r.sum_case);			
			if (r.sum_case>=Seuil_ONLY_CASE) {
				res.significatif=sign_util;
			}
			break;
		case 4:
			res.texte=msprintf("Only one clade");
			break;
			// Manque plein de trucs par rapport à la fonction dans chi2tree...
		default:
			fprintf(stderr, "invalid error %i\n", res.error);
		}
		return res;
	}
	datatype_t p_value;
	int ddl=nb_nodes-1;
	if (r.chi2invalid == 0) {
		if (sign_util) {
			res.significatif=chi2_significatif(ddl, r.chi2);
		}
		p_value=1-gsl_cdf_chisq_P(r.chi2, ddl);
	} else {
		if (texte) {
			res.warning=msprintf("Small sample size correction used");
		}
		// J'ai pas compté dans combien de branches...
		if (ddl == 1) {
			p_value=bilateral(nodes[0].cases,
					  nodes[0].controls,
					  nodes[1].cases,
					  nodes[1].controls);
			if (sign_util) {
				res.significatif=chi2_fisher_significatif(p_value);
			}
		} else {
			p_value = reech_chi2(r.sum_case, r.sum_control,
					     nb_nodes, r.chi2, nodes, th);
			res.warning=msprintfcat(res.warning," (%.6g)",p_value);
			if (sign_util) {
				res.significatif = reech_significatif(p_value);
				if (texte
				    && res.significatif != chi2_significatif(ddl, r.chi2))
				{
					res.warning=msprintfcat(res.warning," Result has changed !");
				}
			}
		}
	}
	res.chi2=r.chi2;
	res.p_val=p_value;
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

datatype_t reech_chi2(int sum_case, int sum_control, int nb_nodes,
		      int chi2_reel, const struct cc *nodes, struct cc *th)
{
	int sum_total=sum_case+sum_control;
	datatype_t p_val=0.0;

	int i,k;
	for(i=0; i<nb_nodes; i++) {
		th[i].controls=((datatype_t)(sum_control*(nodes[i].cases+nodes[i].controls)))/sum_total;
		th[i].cases=((datatype_t)(sum_case*(nodes[i].cases+nodes[i].controls)))/sum_total;
	}

	struct cc clades[nb_nodes];
	for (k=1;k<=PERM; k++){
		random_clades(nb_nodes, nodes, sum_case, sum_control, clades);

		datatype_t chi2=0.0;
		for(i=0; i<nb_nodes; i++) {
			datatype_t t;
			t=clades[i].cases - th[i].cases;
			chi2 += t*t/th[i].cases;
			t=clades[i].controls - th[i].controls;
			chi2 += t*t/th[i].controls;
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

void random_clades(int nb_nodes, const struct cc *nodes,
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
