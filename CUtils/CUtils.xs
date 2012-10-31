#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include "fisher.h"
#include "chisq.h"
#include "double_permutation.h"
#include "resampling.h"
#include "stats.h"

#include "const-c.inc"

MODULE = ALTree::CUtils		PACKAGE = ALTree::CUtils		

INCLUDE: const-xs.inc

############################################################
# fisher.h
############################################################


double
bilateral(a, b, c, d)
	double	a
	double	b
	double	c
	double	d

double
right(a, b, c, d)
	double	a
	double	b
	double	c
	double	d

double
left(a, b, c, d)
	double	a
	double	b
	double	c
	double	d

############################################################
# chisq.h
############################################################

double
pochisq(x, df)
	double	x
	int	df

double
critchi(p, df)
	double	p
	int	df

############################################################
# double_permutation.h
############################################################

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

############################################################
# stats.h
############################################################

void
ClassicalChi2(tabnodes)
	AV * tabnodes
    INIT:
	struct cc *nodes;
	struct classical_chi2_res res;
	int nb_nodes=0;
	int i;
    PPCODE:
	nb_nodes=av_len(tabnodes)+1;
	//fprintf(stderr, "\nSTART(%i, %i)\n", nb_leaves, nb_nodes);
	if (nb_nodes < 1)
	{
	  XSRETURN_UNDEF;
	}
	nodes=(struct cc*)malloc(nb_nodes*sizeof(struct cc));

	for (i=0; i<nb_nodes; i++) {
	  SV* ref=*av_fetch(tabnodes, i, 0);
	  if (!SvROK(ref)) { return XSRETURN_UNDEF; }
	  if (SvTYPE(SvRV(ref)) != SVt_PVHV) { return XSRETURN_UNDEF; }
	  HV* hash=(HV*)SvRV(ref);
	  SV** svp;
	  svp=hv_fetch(hash, "case", 4, 0);
	  if (!svp) { return XSRETURN_UNDEF; }
	  nodes[i].cases=SvNV(*svp);
	  svp=hv_fetch(hash, "control", 7, 0);
	  if (!svp) { return XSRETURN_UNDEF; }
	  nodes[i].controls=SvNV(*svp);
	}

	res=classical_chi2(nb_nodes, nodes);

	EXTEND(SP, 5);
	PUSHs(sv_2mortal(newSVnv(res.chi2)));
	PUSHs(sv_2mortal(newSViv(res.chi2invalid)));
	PUSHs(sv_2mortal(newSViv(res.error)));
	PUSHs(sv_2mortal(newSViv(res.sum_control)));
	PUSHs(sv_2mortal(newSViv(res.sum_case)));

	free(nodes);
