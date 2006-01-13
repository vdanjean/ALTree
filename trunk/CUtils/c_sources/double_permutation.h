//#define CONV "%Le"
//typedef long double datatype_t;

#define CONV "%le"
typedef double datatype_t;

typedef datatype_t *ensemble_t;
typedef ensemble_t *matrice_t;
typedef datatype_t *replicat_t;

ensemble_t alloc_replicat(int nb_chi2);
void free_replicat(ensemble_t rep);

matrice_t alloc_matrice(int nb_sample, int nb_chi2);
void free_matrice(matrice_t mat, int nb_sample, int nb_chi2);

datatype_t calcul(int nb_sample, int nb_chi2, matrice_t mat, replicat_t rep);
