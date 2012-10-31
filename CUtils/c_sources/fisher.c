
#include "fisher.h"
#include <math.h>
#include <stdlib.h>

typedef struct {
	double bilateral;
	double left;
	double right;
} result_t;

static double lngamm(double z)
// Reference: "Lanczos, C. 'A precision approximation
// of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
// Translation of  Alan Miller's FORTRAN-implementation
// See http://lib.stat.cmu.edu/apstat/245
{
  double x = 0;
  x += 0.1659470187408462e-06/(z+7);
  x += 0.9934937113930748e-05/(z+6);
  x -= 0.1385710331296526    /(z+5);
  x += 12.50734324009056     /(z+4);
  x -= 176.6150291498386     /(z+3);
  x += 771.3234287757674     /(z+2);
  x -= 1259.139216722289     /(z+1);
  x += 676.5203681218835     /(z);
  x += 0.9999999999995183;
  return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
}



static double lnfact(double n)
{
  if(n<=1) return(0);
  return(lngamm(n+1));
}

static double lnbico(double n,double k)
{
  return(lnfact(n)-lnfact(k)-lnfact(n-k));
}

static double hyper_323(double n11,double n1_,double n_1,double n)
{
  return(exp(lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1)));
}

static double sn11,sn1_,sn_1,sn,sprob;
static double hyper0(double n11i,double n1_i,double n_1i,double ni)
{
  if(!(n1_i||n_1i||ni))
  {
    if(!(((long long)n11i) % 10 == 0))
    {
      if(n11i==sn11+1)
      {
        sprob *= ((sn1_-sn11)/(n11i))*((sn_1-sn11)/(n11i+sn-sn1_-sn_1));
        sn11 = n11i;
        return sprob;
      }
      if(n11i==sn11-1)
      {
        sprob *= ((sn11)/(sn1_-n11i))*((sn11+sn-sn1_-sn_1)/(sn_1-n11i));
        sn11 = n11i;
        return sprob;
      }
    }
    sn11 = n11i;
  }
  else
  {
    sn11 = n11i;
    sn1_=n1_i;
    sn_1=n_1i;
    sn=ni;
  }
  sprob = hyper_323(sn11,sn1_,sn_1,sn);
  return sprob;
}

static double hyper(double n11)
{
  return(hyper0(n11,0,0,0));
}

static double sleft,sright,sless,slarg;
static double exact(double n11,double n1_,double n_1,double n)
{
  double p,i,j,prob;
  double max=n1_;
  double min;
  if(n_1<max) max=n_1;
  min = n1_+n_1-n;
  if(min<0) min=0;
  if(min==max)
  {
    sless = 1;
    sright= 1;
    sleft = 1;
    slarg = 1;
    return 1;
  }
  prob=hyper0(n11,n1_,n_1,n);
  sleft=0;
  p=hyper(min);
  for(i=min+1; p<0.99999999*prob; i++)
  {
    sleft += p;
    p=hyper(i);
  }
  i--;
  if(p<1.00000001*prob) sleft += p;
  else i--;
  sright=0;
  p=hyper(max);
  for(j=max-1; p<0.99999999*prob; j--)
  {
    sright += p;
    p=hyper(j);
  }
  j++;
  if(p<1.00000001*prob) sright += p;
  else j++;
  if(abs(i-n11)<abs(j-n11))
  {
    sless = sleft;
    slarg = 1 - sleft + prob;
  }
  else
  {
    sless = 1 - sright + prob;
    slarg = sright;
  }
  return prob;
}


/* double left,right,twotail; */
/* double exact22(double n11_,double n12_,double n21_,double n22_) */
/* { */
/*   if(n11_<0) n11_ *= -1; */
/*   if(n12_<0) n12_ *= -1; */
/*   if(n21_<0) n21_ *= -1; */
/*   if(n22_<0) n22_ *= -1; */

/*   double n1_ = n11_+n12_; */
/*   double n_1 = n11_+n21_; */
/*   double n   = n11_ +n12_ +n21_ +n22_; */
/*   double prob=exact(n11_,n1_,n_1,n); */
/*   left    = sless; */
/*   right   = slarg; */
/*   twotail = sleft+sright; */
/*   if(twotail>1) twotail=1; */

/*   "Left   : p-value = "+ left + newline + */
/*   "Right  : p-value = "+ right + newline + */
/*   "2-Tail : p-value = "+ twotail + */
/*   newline +   "------------------------------------------"; */
/* } */

static result_t fisher(double n11, double n12, double n21, double n22)
{
	result_t res;
	double n1_;
	double n_1;
	double n;
	//double prob;

	if(n11<0) n11 *= -1;
	if(n12<0) n12 *= -1;
	if(n21<0) n21 *= -1;
	if(n22<0) n22 *= -1;

	n1_ = n11+n12;
	n_1 = n11+n21;
	n   = n11 +n12 +n21 +n22;
	//prob=
	exact(n11,n1_,n_1,n);
	
	res.bilateral=sleft+sright;
	res.left=sless;
	res.right=slarg;
	if (res.bilateral>1) res.bilateral=1;
	return res;
}

double bilateral(double a, double b, double c, double d)
{
	result_t res=fisher(a,b,c,d);
	//return 1;
	return res.bilateral;
}

double right(double a, double b, double c, double d)
{
	result_t res=fisher(a,b,c,d);
	return res.right;
}

double left(double a, double b, double c, double d)
{
	result_t res=fisher(a,b,c,d);
	return res.left;
}


