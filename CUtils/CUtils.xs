#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include "fisher.h"
#include "chisq.h"
#include "double_permutation.h"

#include "const-c.inc"

MODULE = ALTree::CUtils		PACKAGE = ALTree::CUtils		

INCLUDE: const-xs.inc

double
bilateral(a, b, c, d)
	double	a
	double	b
	double	c
	double	d

double
critchi(p, df)
	double	p
	int	df

double
left(a, b, c, d)
	double	a
	double	b
	double	c
	double	d

double
pochisq(x, df)
	double	x
	int	df

double
right(a, b, c, d)
	double	a
	double	b
	double	c
	double	d

HV *
DoublePermutation(nb_sample, nb_chi2, data)
	int nb_sample
	int nb_chi2
	AV * data
    INIT:
        matrice_t mat;
        replicat_t rep;
        ensemble_t ens;
        int i,j,n;
        datatype_t min;
        AV * ra;
    CODE:
	//fprintf(stderr, "\nSTART(%i, %i)\n", nb_sample, nb_chi2);
        if ((nb_sample <= 0)
	    || (nb_chi2 <= 0)
            || (av_len(data) != nb_sample*nb_chi2-1))
        {
		XSRETURN_UNDEF;
	}
	
        mat=alloc_matrice(nb_sample, nb_chi2);
	rep=alloc_replicat(nb_chi2);
        ens=alloc_ensemble(nb_sample);
	
        n=0;
	for (i=0; i<nb_sample; i++) {
		for (j=0; j<nb_chi2; j++) {
			//fprintf(stderr, "[%i][%i](%i)...\n", i, j, n);
                        /* Attention: on place un réplicat par colonne
                         * (et pas par ligne) */
			mat[j][i]=SvNV(*av_fetch(data, n, 0));
			//fprintf(stderr, "[%i][%i](%i)=%lf\n", i, j, n,mat[j][i]);
			n++;
		}
	}
	min=double_permutation(nb_sample, nb_chi2, mat, rep, ens);

        RETVAL = newHV();
        sv_2mortal((SV *)RETVAL);

	hv_store(RETVAL, "pmin", 4, newSVnv(min), 0);

	ra = (AV *)sv_2mortal((SV *)newAV());
	for (j=0; j<nb_chi2; j++) {
		av_push(ra, newSVnv(rep[j]));
		//fprintf(stderr, "chi2[%i]=%lf\n", j ,rep[j]);
	}	
	hv_store(RETVAL, "chi2", 4, newRV((SV *)ra), 0);

	ra = (AV *)sv_2mortal((SV *)newAV());
	for (j=0; j<nb_sample; j++) {
		av_push(ra, newSVnv(ens[j]));
		//fprintf(stderr, "pmin[%i]=%lf\n", j ,rep[j]);
	}	
	hv_store(RETVAL, "distrib_pmin", 12, newRV((SV *)ra), 0);

        free_matrice(mat, nb_sample, nb_chi2);
        free_replicat(rep);
        free_ensemble(ens);
    OUTPUT:
        RETVAL
