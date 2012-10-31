#include "myrand.h"
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include "debug.h"

#define MT

#ifdef MT
#  include "mt19937ar.h"

static __thread mt19937ar_t rdata;
#  define MYRAND_MAX (0xffffffffUL)
#else
#  define MYRAND_MAX RAND_MAX
#endif


int myrand_init(unsigned long value) {
#ifdef MT
	init_genrand_mt(&rdata, value);
#endif
	return 0;
}

static void init() __attribute__((constructor));
static void init() {
#ifdef MT
	myrand_init(getpid());
#else
	srand(getpid());
#endif
	debug("init rand done");
}


int myrand(int up) {

	for(;;) {
#ifdef MT
		unsigned long r=genrand_int32_mt(&rdata);
#else
		int r=rand();
#endif

		if (r<=MYRAND_MAX-up) {
			return r%up;
		}
		if (r<MYRAND_MAX-(MYRAND_MAX%up)) {
			return r%up;
		}
	}
}
