#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include <c_sources/fisher.h>
#include <c_sources/chisq.h>
#include <c_sources/double_permutation.h>

#include "const-c.inc"

MODULE = Alphy::CUtils		PACKAGE = Alphy::CUtils		

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

SV *
double_permutation(nb_sample, nb_chi2, data)
	int nb_sample
	int nb_chi2
	SV * data
    INIT:
        matrice_t mat;
        replicat_t rep;
        int i,j,n;
        datatype_t min;
	HV * rh;
        AV * ra;
    CODE:
	//fprintf(stderr, "\nSTART(%i, %i)\n", nb_sample, nb_chi2);
        if ((nb_sample <= 0)
	    || (nb_chi2 <= 0)
	    || (!SvROK(data))
            || (SvTYPE(SvRV(data)) != SVt_PVAV)
            || (av_len((AV *)SvRV(data)) != nb_sample*nb_chi2-1))
        {
		XSRETURN_UNDEF;
	}
        rh = (HV *)sv_2mortal((SV *)newHV());
	
        mat=alloc_matrice(nb_sample, nb_chi2);
	rep=alloc_replicat(nb_chi2);
	
        n=0;
	for (i=0; i<nb_sample; i++) {
		for (j=0; j<nb_chi2; j++) {
			//fprintf(stderr, "[%i][%i](%i)...\n", i, j, n);
                        /* Attention: on place un réplicat par colonne
                         * (et pas par ligne) */
			mat[j][i]=SvNV(*av_fetch((AV *)SvRV(data), n, 0));
			//fprintf(stderr, "[%i][%i](%i)=%lf\n", i, j, n,mat[j][i]);
			n++;
		}
	}
	min=calcul(nb_sample, nb_chi2, mat, rep);

	hv_store(rh, "pmin", 4, newSVnv(min), 0);

	ra = (AV *)sv_2mortal((SV *)newAV());
	for (j=0; j<nb_chi2; j++) {
		av_push(ra, newSVnv(rep[j]));
		//fprintf(stderr, "chi2[%i]=%lf\n", j ,rep[j]);
	}
	
	hv_store(rh, "chi2", 4, newRV((SV *)ra), 0);

        free_matrice(mat, nb_sample, nb_chi2);
        free_replicat(rep);

	RETVAL = newRV((SV *)rh);
    OUTPUT:
        RETVAL
