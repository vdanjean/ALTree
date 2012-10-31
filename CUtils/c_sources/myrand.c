#include "myrand.h"
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include "debug.h"

int myrand(int up) {

	static int init=0;
	if (!init) {
		srand(getpid());
		init=1;
		//debug("init rand done");
	}

	for(;;) {
		int r=rand();

		if (r<=RAND_MAX-up) {
			return r%up;
		}
		if (r<RAND_MAX-(RAND_MAX%up)) {
			return r%up;
		}
	}
}
