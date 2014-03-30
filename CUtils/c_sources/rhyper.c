/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2012 The R Core Team
 *  Copyright (C) 2005	The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rhyper(double NR, double NB, double n);
 *
 *  DESCRIPTION
 *
 *    Random variates from the hypergeometric distribution.
 *    Returns the number of white balls drawn when kk balls
 *    are drawn at random from an urn containing nn1 white
 *    and nn2 black balls.
 *
 *  REFERENCE
 *
 *    V. Kachitvichyanukul and B. Schmeiser (1985).
 *    ``Computer generation of hypergeometric random variates,''
 *    Journal of Statistical Computation and Simulation 22, 127-145.
 *
 *    The original algorithm had a bug -- R bug report PR#7314 --
 *    giving numbers slightly too small in case III h2pe
 *    where (m < 100 || ix <= 50) , see below.
 */

#include "mt19937ar.h"
#include "debug.h"
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "rhyper.h"

/* afc(i) :=  ln( i! )	[logarithm of the factorial i.
 *	   If (i > 7), use Stirling's approximation, otherwise use table lookup.
*/

struct afc_data {
    int computed;
    double al[1756];
};

double compute(int n, struct afc_data * __restrict__ data) {
    static long double cur=3628800;
    static int i=11;
    static volatile int mutex=0;

    while (__sync_lock_test_and_set(&mutex, 1)) {
	/* Internal loop with only read to avoid cache line ping-pong
	   on multi-processors */
	while(mutex) {
	    /* spinlock */
	}
     }

     for(; i<=n; i++) {
	cur*=i;
	data->al[i+1]=logl(cur);
     }
     data->computed=n;
     __sync_lock_release(&mutex);
     return data->al[i];
};

static double afc(int i)
{
    double di, value;
    static struct afc_data data = {
	.computed = 10,
	.al = {
	    0.0,
	    0,/*ln(0!)*/
	    0,/*ln(1!)*/
	    0.693147180559945309,/*ln(2!)*/
	    1.791759469228055,/*ln(3!)*/
	    3.17805383034794562,/*ln(4!)*/
	    4.78749174278204599,/*ln(5!)*/
	    6.579251212010101,/*ln(6!)*/
	    8.5251613610654143,/*ln(7!)*/
	    10.6046029027452502,/*ln(8!)*/
	    12.8018274800814696,/*ln(9!)*/
	    15.1044125730755153,/*ln(10!)*/
	}
    };

    if (i < 0) {
      fprintf(stderr, "rhyper.c: afc(i), i=%d < 0 -- SHOULD NOT HAPPEN!\n", i);
      exit(1);
    } else if (i <= data.computed) {
	value = data.al[i + 1];
    } else if (i <= 1754) {
	value = compute(i, &data);
    } else {
	di = i;
	value = (di + 0.5) * log(di) - di + 0.08333333333333 / di
	    - 0.00277777777777 / di / di / di + 0.9189385332;
    }
    return value;
}

#define imin2(a,b) ({ \
	typeof(a) _a=(a); \
	typeof(b) _b=(b); \
	(_a < _b) ? _a : _b ;\
})

#define imax2(a,b) ({ \
	typeof(a) _a=(a); \
	typeof(b) _b=(b); \
	(_a > _b) ? _a : _b ;\
})

#define unif_rand() genrand_real2()

int rhyper(int nn1, int nn2, int kk)
{
    const static double con = 57.56462733;
    const static double deltal = 0.0078;
    const static double deltau = 0.0034;
    const static double scale = 1e25;

    /* extern double afc(int); */

    int i, ix;
    int reject, setup1, setup2;

    double e, f, g, p, r, t, u, v, y;
    double de, dg, dr, ds, dt, gl, gu, nk, nm, ub;
    double xk, xm, xn, y1, ym, yn, yk, alv;

    /* These should become `thread_local globals' : */
    //int ks = -1;
    //int n1s = -1, n2s = -1;

    int k, m;
    int minjx, maxjx, n1, n2;

    double a, d, s, w;
    double tn, xl, xr, kl, kr, lamdl, lamdr, p1, p2, p3;


    /* check parameter validity */

    if (nn1 < 0 || nn2 < 0 || kk < 0 || kk > nn1 + nn2)
	return -1;

    /* if new parameter values, initialize */
    reject = 1;
    //if (nn1 != n1s || nn2 != n2s) {
	setup1 = 1;	setup2 = 1;
    /*} else if (kk != ks) {
	setup1 = 0;	setup2 = 1;
    } else {
	setup1 = 0;	setup2 = 0;
    }*/
    if (setup1) {
	//n1s = nn1;
	//n2s = nn2;
	tn = nn1 + nn2;
	if (nn1 <= nn2) {
	    n1 = nn1;
	    n2 = nn2;
	} else {
	    n1 = nn2;
	    n2 = nn1;
	}
    }
    if (setup2) {
	//ks = kk;
	if (kk + kk >= tn) {
	    k = (int)(tn - kk);
	} else {
	    k = kk;
	}
    }
    if (setup1 || setup2) {
	m = (int) ((k + 1.0) * (n1 + 1.0) / (tn + 2.0));
	minjx = imax2(0, k - n2);
	maxjx = imin2(n1, k);
    }
    /* generate random variate --- Three basic cases */

    if (minjx == maxjx) { /* I: degenerate distribution ---------------- */
	ix = maxjx;
	/* return ix;
	   No, need to unmangle <TSL>*/
	/* return appropriate variate */

	if (kk + kk >= tn) {
	  if (nn1 > nn2) {
	    ix = kk - nn2 + ix;
	  } else {
	    ix = nn1 - ix;
	  }
	} else {
	  if (nn1 > nn2)
	    ix = kk - ix;
	}
	//debug("RHYPER: (%i, %i, %i)=%i", nn1, nn2, kk, ix);
	assert(ix <= nn1);
	assert(kk-ix <= nn2);
	assert(ix <= kk);
	assert(0 <= ix);
	return ix;

    } else if (m - minjx < 10) { /* II: inverse transformation ---------- */
	if (setup1 || setup2) {
	    if (k < n2) {
		w = exp(con + afc(n2) + afc(n1 + n2 - k)
			- afc(n2 - k) - afc(n1 + n2));
	    } else {
		w = exp(con + afc(n1) + afc(k)
			- afc(k - n2) - afc(n1 + n2));
	    }
	}
      L10:
	p = w;
	ix = minjx;
	u = unif_rand() * scale;
      L20:
	if (u > p) {
	    u -= p;
	    p *= (n1 - ix) * (k - ix);
	    ix++;
	    p = p / ix / (n2 - k + ix);
	    if (ix > maxjx)
		goto L10;
	    goto L20;
	}
    } else { /* III : h2pe --------------------------------------------- */

	if (setup1 || setup2) {
	    s = sqrt((tn - k) * k * n1 * n2 / (tn - 1) / tn / tn);

	    /* remark: d is defined in reference without int. */
	    /* the truncation centers the cell boundaries at 0.5 */

	    d = (int) (1.5 * s) + .5;
	    xl = m - d + .5;
	    xr = m + d + .5;
	    a = afc(m) + afc(n1 - m) + afc(k - m) + afc(n2 - k + m);
	    kl = exp(a - afc((int) (xl)) - afc((int) (n1 - xl))
		     - afc((int) (k - xl))
		     - afc((int) (n2 - k + xl)));
	    kr = exp(a - afc((int) (xr - 1))
		     - afc((int) (n1 - xr + 1))
		     - afc((int) (k - xr + 1))
		     - afc((int) (n2 - k + xr - 1)));
	    lamdl = -log(xl * (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1));
	    lamdr = -log((n1 - xr + 1) * (k - xr + 1) / xr / (n2 - k + xr));
	    p1 = d + d;
	    p2 = p1 + kl / lamdl;
	    p3 = p2 + kr / lamdr;
	}
      L30:
	u = unif_rand() * p3;
	v = unif_rand();
	if (u < p1) {		/* rectangular region */
	    ix = (int) (xl + u);
	} else if (u <= p2) {	/* left tail */
	    ix = (int) (xl + log(v) / lamdl);
	    if (ix < minjx)
		goto L30;
	    v = v * (u - p1) * lamdl;
	} else {		/* right tail */
	    ix = (int) (xr - log(v) / lamdr);
	    if (ix > maxjx)
		goto L30;
	    v = v * (u - p2) * lamdr;
	}

	/* acceptance/rejection test */

	if (m < 100 || ix <= 50) {
	    /* explicit evaluation */
	    /* The original algorithm (and TOMS 668) have
		   f = f * i * (n2 - k + i) / (n1 - i) / (k - i);
	       in the (m > ix) case, but the definition of the
	       recurrence relation on p134 shows that the +1 is
	       needed. */
	    f = 1.0;
	    if (m < ix) {
		for (i = m + 1; i <= ix; i++)
		    f = f * (n1 - i + 1) * (k - i + 1) / (n2 - k + i) / i;
	    } else if (m > ix) {
		for (i = ix + 1; i <= m; i++)
		    f = f * i * (n2 - k + i) / (n1 - i + 1) / (k - i + 1);
	    }
	    if (v <= f) {
		reject = 0;
	    }
	} else {
	    /* squeeze using upper and lower bounds */
	    y = ix;
	    y1 = y + 1.0;
	    ym = y - m;
	    yn = n1 - y + 1.0;
	    yk = k - y + 1.0;
	    nk = n2 - k + y1;
	    r = -ym / y1;
	    s = ym / yn;
	    t = ym / yk;
	    e = -ym / nk;
	    g = yn * yk / (y1 * nk) - 1.0;
	    dg = 1.0;
	    if (g < 0.0)
		dg = 1.0 + g;
	    gu = g * (1.0 + g * (-0.5 + g / 3.0));
	    gl = gu - .25 * (g * g * g * g) / dg;
	    xm = m + 0.5;
	    xn = n1 - m + 0.5;
	    xk = k - m + 0.5;
	    nm = n2 - k + xm;
	    ub = y * gu - m * gl + deltau
		+ xm * r * (1. + r * (-0.5 + r / 3.0))
		+ xn * s * (1. + s * (-0.5 + s / 3.0))
		+ xk * t * (1. + t * (-0.5 + t / 3.0))
		+ nm * e * (1. + e * (-0.5 + e / 3.0));
	    /* test against upper bound */
	    alv = log(v);
	    if (alv > ub) {
		reject = 1;
	    } else {
				/* test against lower bound */
		dr = xm * (r * r * r * r);
		if (r < 0.0)
		    dr /= (1.0 + r);
		ds = xn * (s * s * s * s);
		if (s < 0.0)
		    ds /= (1.0 + s);
		dt = xk * (t * t * t * t);
		if (t < 0.0)
		    dt /= (1.0 + t);
		de = nm * (e * e * e * e);
		if (e < 0.0)
		    de /= (1.0 + e);
		if (alv < ub - 0.25 * (dr + ds + dt + de)
		    + (y + m) * (gl - gu) - deltal) {
		    reject = 0;
		}
		else {
		    /* * Stirling's formula to machine accuracy
		     */
		    if (alv <= (a - afc(ix) - afc(n1 - ix)
				- afc(k - ix) - afc(n2 - k + ix))) {
			reject = 0;
		    } else {
			reject = 1;
		    }
		}
	    }
	}
	if (reject)
	    goto L30;
    }

    /* return appropriate variate */

    if (kk + kk >= tn) {
	if (nn1 > nn2) {
	    ix = kk - nn2 + ix;
	} else {
	    ix = nn1 - ix;
	}
    } else {
	if (nn1 > nn2)
	    ix = kk - ix;
    }
    //debug("RHYPER: (%i, %i, %i)=%i", nn1, nn2, kk, ix);
    assert(ix <= nn1);
    assert(kk-ix <= nn2);
    assert(ix <= kk);
    assert(0 <= ix);
    return ix;
}

#if TEST_AFC
static double origafc(int i)
{
    const static double al[9] =
    {
	0.0,
	0.0,/*ln(0!)=ln(1)*/
	0.0,/*ln(1!)=ln(1)*/
	0.69314718055994530941723212145817,/*ln(2) */
	1.79175946922805500081247735838070,/*ln(6) */
	3.17805383034794561964694160129705,/*ln(24)*/
	4.78749174278204599424770093452324,
	6.57925121201010099506017829290394,
	8.52516136106541430016553103634712
	/*, 10.60460290274525022841722740072165*/
    };
    double di, value;

    if (i < 0) {
      fprintf(stderr, "rhyper.c: afc(i), i=%d < 0 -- SHOULD NOT HAPPEN!\n", i);
      exit(1);
    } else if (i <= 7) {
	value = al[i + 1];
    } else {
	di = i;
	value = (di + 0.5) * log(di) - di + 0.08333333333333 / di
	    - 0.00277777777777 / di / di / di + 0.9189385332;
    }
    return value;
}

static double afc2(int n)
{
	static const double logpi=__builtin_log(M_PI);
	
	return n*log(n)-n+log(n*(1+4*n*(1+2*n)))/6+logpi/2;
}

static double afc3(int n)
{
	static const double logpi=__builtin_log(M_PI);
	
	return n*log(n)-n+log(1+1/(2*n)+1/(8*n*n))/6+log(2*n)/2+logpi/2;
}

static double afc4(int n)
{
	static const double logpi=__builtin_log(M_PI);
	static const double log2=__builtin_log(2);
	double logn=log(n);
	
	return n*logn-n+log(1+1/(2*n)+1/(8*n*n))/6+(logn+(logpi+log2))/2;
}

static double afc5(int n)
{
	static long long int cur=1;
	static int i=1;

	for(; i<=n; i++) {
		cur*=i;
	}
	//printf(" %lli %i %i ", cur, i, n);
	return log(cur);
}

static double afc6(int n)
{
	static long double cur=1;
	static int i=1;

	for(; i<=n; i++) {
		cur*=i;
	}
	//printf(" %lli %i %i ", cur, i, n);
	return logl(cur);
}

static double afc7(int n)
{
	static long double cur=1;
	static int i=1;

	for(; i<=n; i++) {
		cur*=i;
	}
	//printf(" %lli %i %i ", cur, i, n);
	printf("\t%.18Lg, /* ln(%i!) = ln(%.0Lf) */\n",logl(cur),n,cur);
	return logl(cur);
}

static void compare(int k) {
	int i;
	printf("           %20s / %20s / %23s / %23s / %23s / %23sg / \n",
		"ref=exact(long double)", "my", "orig-ref", "orig-my", "my-ref", "exact(double)-ref");
	for (i=1; i<=k; i++) {
		double ref=afc6(i);
		double ref2=afc(i);
		printf("log %4i! = %20.17lg / %20.17lg / %13.7lg / %13.7lg / %13.7lg / %13.7lg / %13.7lg / %13.7lg / %13.7lg\n",
			i, ref, ref2, origafc(i)-ref, origafc(i)-ref2, afc(i)-ref, afc5(i)-ref, afc2(i)-ref, afc3(i)-ref, afc4(i)-ref);
	}
	for (i=1; i<=50; i++) {
		afc7(i);
	}
	return;
}

int main(int argc, char**argv) {

	compare(1755);
	return 0;
}
#endif
