#ifndef _RESAMPLING_H
#define _RESAMPLING_H

#include "datatype.h"
#include "stats.h"

struct tree {
	int nb_leaves;
	int nb_nodes;
	int max_depth;
	int *ld;
	int *lp;
	int *np;
};

int resampling_chi2(const struct tree *tree, const struct cc *lcc,
		    int prolonge, int nb_permutations,
		    datatype_t *results, int parallel);

#endif
