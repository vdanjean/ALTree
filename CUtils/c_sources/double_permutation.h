#ifndef _DOUBLE_PERMUTATION_H
#define _DOUBLE_PERMUTATION_H

#include "datatype.h"

typedef datatype_t *ensemble_t;
typedef ensemble_t *matrice_t;
typedef datatype_t *replicat_t;

ensemble_t alloc_replicat(int nb_chi2);
void free_replicat(ensemble_t rep);

ensemble_t alloc_ensemble(int nb_sample);
void free_ensemble(ensemble_t ens);

matrice_t alloc_matrice(int nb_sample, int nb_chi2);
void free_matrice(matrice_t mat, int nb_sample, int nb_chi2);

datatype_t calcul(int nb_sample, int nb_chi2, matrice_t mat, replicat_t rep);
datatype_t double_permutation(int nb_sample, int nb_chi2, matrice_t mat,
			      replicat_t rep, ensemble_t ens_min_pval);

#endif

