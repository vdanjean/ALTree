
#include "debug.h"
#include "resampling.h"
#include <strings.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/types.h>
#include <unistd.h>
#include "myrand.h"

struct chi2s {
	int nb_leaves;
	struct cc*leaves;
	struct cc*th;
	int nb_nodes;
	int *stashed_nodes;
	int *next_nodes;
	struct cc*nodes;
};

static void compute_chi2s(const struct tree *tree, const struct cc *lcc,
			  struct chi2s *temp,
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
		struct calcul_chi2_res r=calcul_chi2(nb_leaves, temp->leaves,
						     0, 0, temp->th);
		assert(r.error==0);
		results[depth-1]=r.chi2;
	}
}

struct paractl {
	volatile int permutation;
	int nb_permutations;
	const struct tree *tree;
	const struct cc *lcc;
	int cases;
	int controls;
	int prolonge;
	datatype_t *results;
};

struct parainfo {
	struct paractl *ctl;
	int permutation;
};

struct memory {
	struct cc *rand_lcc;
	struct chi2s temp;
};

static struct memory *mem_alloc(const struct tree *tree)
{
	struct memory *m=malloc(sizeof(struct memory));

	m->rand_lcc=malloc(tree->nb_leaves*sizeof(struct cc));
	m->temp=(struct chi2s){
		.nb_leaves=tree->nb_leaves,
		.leaves=malloc(tree->nb_leaves*sizeof(struct cc)),
		.th=malloc(tree->nb_leaves*sizeof(struct cc)),
		.nb_nodes=tree->nb_nodes,
		.stashed_nodes=malloc(tree->nb_nodes*sizeof(int)),
		.next_nodes=malloc(tree->nb_nodes*sizeof(int)),
		.nodes=malloc(tree->nb_nodes*sizeof(struct cc)),
	};

	return m;
}

static void free_alloc(struct memory *m)
{
	free(m->temp.leaves);
	free(m->temp.th);
	free(m->temp.stashed_nodes);
	free(m->temp.next_nodes);
	free(m->temp.nodes);
	free(m->rand_lcc);
}

static void *resampling_worker(void* param) {
	struct parainfo *p=(struct parainfo *)param;
	
	/* Copy values for ro access, to avoid ping-pong with cache lines
	   when updating permutation */
	struct paractl data= *(p->ctl);
	int permutation=p->permutation;

	struct memory *mem=mem_alloc(data.tree);

	myrand_init(pthread_self()+getpid());

	while (permutation < data.nb_permutations) {
		//debug("thread %i handle perm %i", p->permutation, permutation);
		random_clades(data.tree->nb_leaves, data.lcc,
			      data.cases, data.controls, mem->rand_lcc);
		compute_chi2s(data.tree, mem->rand_lcc, &mem->temp, data.prolonge,
			      data.results+(permutation * data.tree->max_depth));
		permutation=__sync_fetch_and_add(&p->ctl->permutation, 1);
	}

	free_alloc(mem);

	return NULL;
}

/* From http://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine */
static int nbproc()
{
	int numCPU=0;
#ifdef __linux
	numCPU = sysconf( _SC_NPROCESSORS_ONLN );
#elif defined(__bsd)
	int mib[4];
	size_t len = sizeof(numCPU); 

	/* set the mib for hw.ncpu */
	mib[0] = CTL_HW;
	mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;

	/* get the number of CPUs from the system */
	sysctl(mib, 2, &numCPU, &len, NULL, 0);

	if( numCPU < 1 ) {
	     mib[1] = HW_NCPU;
	     sysctl( mib, 2, &numCPU, &len, NULL, 0 );

	     if( numCPU < 1 ) {
	          numCPU = 1;
	     }
	}
#else
#  warning no support on this plate-form. Patch welcome.
#endif
	return numCPU;
}

int resampling_chi2(const struct tree *tree, const struct cc *lcc, int prolonge,
		    int nb_permutations, datatype_t *results, int parallel)
{
#if defined(DEBUG) && 0 
	FILE* dump=fopen("/tmp/dump", "w+");
	fwrite(tree, sizeof(struct tree), 1, dump);
	fwrite(tree->ld, sizeof(int), tree->nb_leaves, dump);
	fwrite(tree->lp, sizeof(int), tree->nb_leaves, dump);
	fwrite(tree->np, sizeof(int), tree->nb_nodes, dump);
	fwrite(lcc, sizeof(struct cc), tree->nb_leaves, dump);
	fwrite(&prolonge, sizeof(int), 1, dump);
	fwrite(&nb_permutations, sizeof(int), 1, dump);
	fclose(dump);
#endif
	char* envvar=getenv("ALTREE_PARALLEL");
	if (envvar) {
		parallel=atoi(envvar);
	}
	if (parallel == -1) {
		parallel=nbproc();
	}
	if (parallel < 0) {
		parallel=0;
	}
	debug("parallel=%i", parallel);

	struct memory *mem=mem_alloc(tree);
	compute_chi2s(tree, lcc, &mem->temp, prolonge, results);
	int i;
	int cases=0;
	int controls=0;
	for(i=0; i<tree->nb_leaves; i++) {
		cases += lcc[i].cases;
		controls += lcc[i].controls;
	}
	if (!parallel) {
		for(i=0; i<nb_permutations; i++) {
			results += tree->max_depth;
			random_clades(tree->nb_leaves, lcc, cases, controls, mem->rand_lcc);
			compute_chi2s(tree, mem->rand_lcc, &mem->temp, prolonge, results);
		}
	} else {
		struct paractl ctl = {
			.permutation=parallel,
			.nb_permutations=nb_permutations,
			.tree=tree,
			.lcc=lcc,
			.cases=cases,
			.controls=controls,
			.prolonge=prolonge,
			.results=results+tree->max_depth,
		};

		struct parainfo infos[parallel];
		pthread_t tids[parallel];

		for(i=0; i<parallel; i++) {
			infos[i].ctl=&ctl;
			infos[i].permutation=i;
			pthread_create(&tids[i], NULL, &resampling_worker, &infos[i]);
		}
		for(i=0; i<parallel; i++) {
			pthread_join(tids[i], NULL);
		}
	}
	free_alloc(mem);

	return 0;
}

#ifdef MAIN_PROG
int main(int argc, char *argv[])
{
	int resampling_chi2(const struct tree *tree, const struct cc *lcc, int prolonge,
			    int nb_permutations, datatype_t *results, int parallel);
	struct tree tree;

	FILE* dump=fopen("/tmp/dump.read", "r");
	fread(&tree, sizeof(struct tree), 1, dump);

	int ld[tree.nb_leaves];
	int lp[tree.nb_leaves];
	int np[tree.nb_nodes];
	struct cc lcc[tree.nb_leaves];
	int prolonge;
	int nb_permutations;

	tree.ld=ld;
	tree.lp=lp;
	tree.np=np;

	fread(tree.ld, sizeof(int), tree.nb_leaves, dump);
	fread(tree.lp, sizeof(int), tree.nb_leaves, dump);
	fread(tree.np, sizeof(int), tree.nb_nodes, dump);
	fread(lcc, sizeof(struct cc), tree.nb_leaves, dump);
	fread(&prolonge, sizeof(int), 1, dump);
	fread(&nb_permutations, sizeof(int), 1, dump);
	nb_permutations=10;
	fclose(dump);

	datatype_t results[(nb_permutations+1)*tree.max_depth];
	//bzero(results, (nb_permutations+1)*tree.max_depth*sizeof(datatype_t));

	resampling_chi2(&tree, lcc, prolonge,
			nb_permutations, results, 0);

	int i,j;
	datatype_t *r=results;
	for(i=0; i<=nb_permutations; i++) {
		for(j=0; j<tree.max_depth; j++) {
			printf("\t"CONV, *(r++));
		}
		printf("\n");
	}

	/* ensemble_t ens_min_pval; */

	/* ens_min_pval=alloc_ensemble(nb_sample); */
	
	/* min=double_permutation(nb_sample, nb_chi2, mat, rep, ens_min_pval); */

	/* free_ensemble(ens_min_pval); */

	/* for (j=0; j<nb_chi2; j++) { */
	/* 	printf("chi2 niveau %i, pval nc "CONV"\n", j+1, rep[j]); */
	/* } */
	/* for (i=0; i<nb_sample; i++) { */
	/* 	printf("sample %i, pval min "CONV"\n", i, ens_min_pval[i]); */
	/* } */
	/* printf("pmin corrigé: "CONV"\n", min); */

	/* free_matrice(mat, nb_sample, nb_chi2); */
	/* free_replicat(rep); */

	return 0;
}
#endif
