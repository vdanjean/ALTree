
#include "debug.h"
#include "resampling.h"
#include <strings.h>


struct chi2s {
	int nb_leaves;
	struct cc*leaves;
	int nb_nodes;
	int *stashed_nodes;
	int *next_nodes;
	struct cc*nodes;
};

static void compute_chi2s(struct tree *tree, struct cc *lcc, struct chi2s *temp,
			  int prolonge, datatype_t *results)
{
	int first_leaf=0;

	int nb_next_nodes=0;
	int nb_stashed_nodes;
	
	assert(temp->nb_leaves==tree->nb_leaves);
	assert(temp->nb_nodes==tree->nb_nodes);

	bzero(temp->nodes, tree->nb_nodes*sizeof(struct cc));
	bzero(temp->stashed_nodes, tree->nb_nodes*sizeof(int));

	//debug("max_depth=%i, nb_leaves=%i, nb_nodes=%i",
	//	tree->max_depth, tree->nb_leaves, tree->nb_nodes);

	int depth;
	for (depth=tree->max_depth; depth>0; depth--){
		struct cc *leaves=temp->leaves;
		int nb_leaves=0;
		int next_first_leaf;

		int i;
		nb_stashed_nodes=0;
		for(i=0; i<nb_next_nodes; i++) {
			int id=temp->next_nodes[i];
			temp->leaves[nb_leaves++]=temp->nodes[id];
			int parent=tree->np[id];
			if (parent != -1) {
				if (!temp->stashed_nodes[parent]) {
					temp->stashed_nodes[parent]=1;
					temp->next_nodes[nb_stashed_nodes++]=parent;
				}
				temp->nodes[parent].cases+=temp->nodes[id].cases;
				temp->nodes[parent].controls+=temp->nodes[id].controls;
			}
		}

		for(i=first_leaf; i<tree->nb_leaves && tree->ld[i] == depth; i++) {
			temp->leaves[nb_leaves++]=lcc[i];
			int parent=tree->lp[i];
			if (parent != -1) {
				if (!temp->stashed_nodes[parent]) {
					temp->stashed_nodes[parent]=1;
					temp->next_nodes[nb_stashed_nodes++]=parent;
				}
				temp->nodes[parent].cases+=lcc[i].cases;
				temp->nodes[parent].controls+=lcc[i].controls;
			}
		}
		assert(nb_stashed_nodes <= tree->nb_leaves);
		next_first_leaf=i;
		if (prolonge == 1) {
			for(; i<tree->nb_leaves; i++) {
				temp->leaves[nb_leaves++]=lcc[i];
			}
		}
		nb_next_nodes=nb_stashed_nodes;
		first_leaf=next_first_leaf;

		//debug("depth=%i, ddl=%i", depth, nb_leaves-1);
		struct calcul_chi2_res r=calcul_chi2(nb_leaves, leaves, 0, 0);
		assert(r.error==0);
		results[depth-1]=r.chi2;
	}
}

int resampling_chi2(struct tree *tree, struct cc *lcc, int prolonge,
		    int nb_permutations, datatype_t *results)
{

	struct cc rand_lcc[tree->nb_leaves];

	struct cc leaves[tree->nb_leaves];
	struct cc nodes[tree->nb_nodes];
	int stashed_nodes[tree->nb_nodes];
	int next_nodes[tree->nb_nodes];
	struct chi2s temp={
		.nb_leaves=tree->nb_leaves,
		.leaves=leaves,
		.nb_nodes=tree->nb_nodes,
		.stashed_nodes=stashed_nodes,
		.next_nodes=next_nodes,
		.nodes=nodes,
	};
	       
	compute_chi2s(tree, lcc, &temp, prolonge, results);
	int i;
	int cases=0;
	int controls=0;
	for(i=0; i<tree->nb_leaves; i++) {
		cases += lcc[i].cases;
		controls += lcc[i].controls;
	}
	for(i=0; i<nb_permutations; i++) {
		results += tree->max_depth;
		random_clades(tree->nb_leaves, lcc, cases, controls, rand_lcc);
		compute_chi2s(tree, rand_lcc, &temp, prolonge, results);
	}

	return 0;
}
