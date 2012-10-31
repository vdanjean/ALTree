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
                        /* Attention: on place un r�plicat par colonne
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

void
CalculChi2(tabnodes, ddl, test_results, sign_util)
	AV * tabnodes
        int ddl
        HV * test_results
        int sign_util
    INIT:
	struct cc *nodes;
	int nb_nodes;
	struct calcul_chi2_res res;
	int i;
    PPCODE:
	nb_nodes=ddl+1;
	if (av_len(tabnodes) != ddl)
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

	res=calcul_chi2(nb_nodes, nodes, sign_util, 1);

	if (res.texte) {
	  hv_store(test_results, "texte", 5,
		   newSVpv(res.texte, 0),0);
	  free(res.texte);
	}
	if (sign_util) {
	  hv_store(test_results, "sign", 4,
		   newSViv(res.significatif),0);
	}
	if (res.warning) {
	  hv_store(test_results, "warning", 7,
		   newSVpv(res.warning, 0),0);
	  free(res.warning);
	}
	if (res.error==0) {
	  hv_store(test_results, "chi2", 4,
		   newSVnv(res.chi2),0);
	  hv_store(test_results, "p_val", 5,
		   newSVnv(res.p_val),0);
	}

	free(nodes);

	EXTEND(SP, 2);
	PUSHs(sv_2mortal(newSVnv(res.p_val)));
	PUSHs(sv_2mortal(newSViv(res.significatif)));

void
DefinitionPChi2(p, pprop)
	SV* p
	SV* pprop
    INIT:
        double cp;
	double cpprop;
    PPCODE:
	if (!SvOK(p) || !SvNOK(p)) {
	  cp=-1;
	} else {
	  cp=SvNV(p);
	}
	if (!SvOK(pprop) || !SvNOK(pprop)) {
	  cpprop=-1;
	} else {
	  cpprop=SvNV(pprop);
	}
	definition_p_chi2(cp, cpprop);

int
Chi2Significatif(ddl, chi2)
        int ddl
        double chi2
    CODE:
        RETVAL=chi2_significatif(ddl, chi2);
    OUTPUT:
        RETVAL

int
Chi2FisherSignificatif(pvalue)
        double pvalue
    CODE:
        RETVAL=chi2_fisher_signigicatif(pvalue);
    OUTPUT:
        RETVAL

double
ReechChi2(sum_case, sum_control, nb_nodes, chi2_reel, clades)
        int sum_case
        int sum_control
        int nb_nodes
        double chi2_reel
        AV* clades
    INIT:
	struct cc *nodes;
	int case_avail;
	int i;
    CODE:
	if (av_len(clades)+1 != nb_nodes)
	{
	  XSRETURN_UNDEF;
	}
	nodes=(struct cc*)malloc(nb_nodes*sizeof(struct cc));

	case_avail=sum_case;
	for (i=0; i<nb_nodes; i++) {
	  SV* val=*av_fetch(clades, i, 0);
	  if (!SvIOK(val)) { return XSRETURN_UNDEF; }
	  int nb=SvIV(val);
	  if (nb <= case_avail) {
	    case_avail-=nb;
	    nodes[i].cases=nb;
	    nodes[i].controls=0;
	  } else if (case_avail == 0) {
	    nodes[i].cases=0;
	    nodes[i].controls=nb;
	  } else {
	    nodes[i].cases=case_avail;
	    nodes[i].controls=nb - case_avail;
	    case_avail=0;
	  }
	}

	RETVAL=reech_chi2(sum_case, sum_control, nb_nodes, chi2_reel, nodes);

	free(nodes);
    OUTPUT:
        RETVAL

int
ReechSignificatif(p_val)
        double p_val
    CODE:
        RETVAL=reech_significatif(p_val);
    OUTPUT:
        RETVAL
