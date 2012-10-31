#ifndef _STATS_H
#define _STATS_H

#include "datatype.h"

#define SAMPLESIZE 5

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

struct classical_chi2_res classical_chi2(int nb_nodes, struct cc *nodes);

#endif

