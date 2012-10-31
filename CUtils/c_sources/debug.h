#ifndef _DEBUG_H
#define _DEBUG_H

//#define DEBUG

#ifdef DEBUG
#  include <stdio.h>
#  define debug(str,...) fprintf(stderr, str "\n", ##__VA_ARGS__)
#else
#  define debug(str,...) ((void)0)
#  define NDEBUG
#endif
#include <assert.h>

#endif
