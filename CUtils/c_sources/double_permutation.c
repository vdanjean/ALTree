
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "double_permutation.h"

#define CALC_PVAL(count, nb_sample) \
     (((datatype_t)(count-1))/* On s'enlève soi-même*/ \
      /(nb_sample))


int read_matrice(matrice_t mat, int nb_sample, int nb_chi2)
{
	int i, j;
	datatype_t d;
	int res;
	for (i=0; i<nb_sample; i++) {
		for (j=0; j<nb_chi2; j++) {
			res=scanf(CONV, &d);
			if (res!=1) {
				fprintf(stderr, "Erreur de lecture. Probablement pas assez de données\n");
				exit(1);
			}
			/* Attention: on place un réplicat par colonne
			 * (et pas par ligne) */
			mat[j][i]=d;
		}
	}
	return 0;
}

ensemble_t alloc_replicat(int nb_chi2)
{
	ensemble_t rep=NULL;	
	rep=malloc(nb_chi2*sizeof(datatype_t));
	if (rep==NULL) {
		goto err;
	}
	return rep;
 err:
	fprintf(stderr, "Erreur d'allocation mémoire. Aborting\n");
	exit(1);
}

void free_replicat(ensemble_t rep)
{
	free(rep);
}

ensemble_t alloc_ensemble(int nb_sample)
{
	ensemble_t ens=NULL;	
	ens=malloc(nb_sample*sizeof(datatype_t));
	if (ens==NULL) {
		goto err;
	}
	return ens;
 err:
	fprintf(stderr, "Erreur d'allocation mémoire. Aborting\n");
	exit(1);
}

void free_ensemble(ensemble_t ens)
{
	free(ens);
}

matrice_t alloc_matrice(int nb_sample, int nb_chi2)
{
	int i;
	matrice_t mat=NULL;
	mat=malloc(nb_chi2*sizeof(ensemble_t));
	if (mat==NULL)
		goto err;
	for (i=0; i<nb_chi2; i++) {
		mat[i]=alloc_ensemble(nb_sample);
	}

	return mat;
 err:
	fprintf(stderr, "Erreur d'allocation mémoire. Aborting\n");
	exit(1);

}

void free_matrice(matrice_t mat, int nb_sample, int nb_chi2)
{
	int i;
	for (i=0; i<nb_chi2; i++) {
		free_ensemble(mat[i]);
	}
	free(mat);
}

inline static int count_superieur(ensemble_t ens, datatype_t val_ref,
				  int nb_sample)
{
	int i, count=0;
	for (i=0; i< nb_sample; i++) {
		if (ens[i]>=val_ref) {
			count++;
		}
	}
	//printf( "count=%i\n", count);
	return count;
}

inline static int count_inferieur(ensemble_t ens, datatype_t val_ref,
				  int nb_sample)
{
	int i, count=0;
	for (i=0; i< nb_sample; i++) {
		if (ens[i]<=val_ref) {
			count++;
		}
	}
	//printf( "count=%i\n", count);
	return count;
}

inline static datatype_t pval_min(replicat_t rep, int nb_chi2)
{
	int i;
	datatype_t ret=rep[0];
	for (i=1; i<nb_chi2; i++) {
		if (rep[i]<ret) {
			ret=rep[i];
		}
	}
	return ret;
}

datatype_t calcul(int nb_sample, int nb_chi2, matrice_t mat, replicat_t rep)
{
	datatype_t min;
	ensemble_t ens_min_pval;

	ens_min_pval=alloc_ensemble(nb_sample);
	
	min=double_permutation(nb_sample, nb_chi2, mat, rep, ens_min_pval);

	free_ensemble(ens_min_pval);
	return min;
}

datatype_t double_permutation(int nb_sample, int nb_chi2, matrice_t mat,
			      replicat_t rep, ensemble_t ens_min_pval)
{
	int i, j;
	datatype_t min;
	datatype_t local[nb_chi2];
	
	i=0;
	for (j=0; j<nb_chi2; j++) {
		rep[j]=CALC_PVAL(count_superieur(mat[j], mat[j][i], nb_sample),
				 nb_sample);
		//printf("cal rep[%i]=%lf\n", j, rep[j]);
	}
	/* i is still 0 here */
	ens_min_pval[i]=pval_min(rep, nb_chi2);
	//printf("pmin for sample %i: "CONV"\n", i, ens_min_pval[i]);

	for (i=1; i<nb_sample; i++) {
		for (j=0; j<nb_chi2; j++) {
			local[j]=CALC_PVAL(count_superieur(mat[j], mat[j][i], 
							   nb_sample),
					   nb_sample);
		}
		ens_min_pval[i]=pval_min(local, nb_chi2);
		//printf("pmin for sample %i: "CONV"\n", i, ens_min_pval[i]);
	}
	min=CALC_PVAL(count_inferieur(ens_min_pval, ens_min_pval[0], 
				      nb_sample),
		      nb_sample);
	return min;
}


#ifdef MAIN_PROG
int main(int argc, char *argv[])
{
	int nb_sample, nb_chi2;
	matrice_t mat;
	replicat_t rep;
	int j,i;
	datatype_t min;

	nb_sample=atoi(argv[1]);
	nb_chi2=atoi(argv[2]);
	
	mat=alloc_matrice(nb_sample, nb_chi2);
	rep=alloc_replicat(nb_chi2);

	read_matrice(mat, nb_sample, nb_chi2);
	//printf("Matrice lue\n");
	for(i=0; i<nb_sample; i++) {
		for(j=0; j<nb_chi2; j++) {
			printf("\t"CONV, mat[j][i]);
		}
		printf("\n");
	}

	ensemble_t ens_min_pval;

	ens_min_pval=alloc_ensemble(nb_sample);
	
	min=double_permutation(nb_sample, nb_chi2, mat, rep, ens_min_pval);

	free_ensemble(ens_min_pval);

	for (j=0; j<nb_chi2; j++) {
		printf("chi2 niveau %i, pval nc "CONV"\n", j+1, rep[j]);
	}
	for (i=0; i<nb_sample; i++) {
		printf("sample %i, pval min "CONV"\n", i, ens_min_pval[i]);
	}
	printf("pmin corrigé: "CONV"\n", min);

	free_matrice(mat, nb_sample, nb_chi2);
	free_replicat(rep);

	return 0;
}
#endif
