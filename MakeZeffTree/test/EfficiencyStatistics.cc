#include <math.h>
#include "EfficiencyStatistics.h"
#include "TRandom.h"
#include "TF1.h"

double EfficiencyStatFunc(double* x, double* param); // two params!

/// SAFE FOR ZERO DENOM????

static const int CrossOver = 20;
static double factorial(int n) { 
  double val=1.0;
  for (int i=1; i<=n; i++) val*=i;
  return val;
}

static double logfactorial(int n) {
  const double ln2pi=0.5*log(2*M_PI);
  if (n==0) return 0;
  if (n<15) return log(factorial(n));
  double lf= ln2pi+(n+0.5)*log(n) - n + 1.0/(12.0*n) - 1.0/(360.0*n*n*n);
  //printf("%f %f\n",lf,log(TMath::Factorial(n)));
  return lf;
}    
static double prefactorCalc(int num, int den) {
  if (den==0) {
    return  0;
  } else if (den<CrossOver) {
    return factorial(den+1)/(factorial(num)*factorial(den-num));
  } else {
    return logfactorial(den+1)-logfactorial(num)-logfactorial(den-num);
  } 
}



EfficiencyStatistics::EfficiencyStatistics(double eff0, int denom) {
  eff0_=eff0;
  if (denom<=0) denom_=0;
  else denom_=denom;
  numer_=int(eff0*denom_+0.5);
  prefactor_=prefactorCalc(numer_,denom_);
  rms_=sqrt(eff0_*(1-eff0_))/sqrt(std::max(1,denom_));
  random_=gRandom;
}

EfficiencyStatistics::EfficiencyStatistics(int numer, int denom) {
  if (denom<=0) eff0_=0;
  else eff0_=numer*1.0/denom;
  if (denom<=0) denom_=0;
  else denom_=denom;
  numer_=numer;
  rms_=sqrt(eff0_*(1-eff0_))/sqrt(std::max(1,denom_));
  prefactor_=prefactorCalc(numer_,denom_);
  random_=gRandom;
}

void EfficiencyStatistics::setRandom(TRandom* rand) {
  random_=rand;
}

double EfficiencyStatistics::prob(double x) const {
  if (denom_==0) { return 0; }
  else if (denom_<CrossOver) {
    return prefactor_*pow(x,numer_)*pow(1-x,denom_-numer_);
  } else {
    double b=prefactor_;
    if (x>0.0) b+=numer_*log(x);
    if (x<1.0) b+=(denom_-numer_)*log(1-x);
    return exp(b);
  }
}

double EfficiencyStatFunc(double* x, double* params) {
  static int numer=0;
  static int denom=0;
  static double pref=0;

  if (int(params[0])!=numer || int(params[1])!=denom) {
    numer=int(params[0]);
    denom=int(params[1]);
    pref=prefactorCalc(numer,denom);
    //    printf("n=%d d=%d\n",numer,denom);
  }
  if (denom==0) return 0;
  else if (denom<CrossOver) {
    return pref*pow(*x,numer)*pow(1-*x,denom-numer);
  } if (denom<120) {
    double b=exp(pref);
    if (*x>0.0) b*=pow(*x,numer);
    if (*x<1.0) b*=pow(1-*x,denom-numer);
    return b;
  } else {
    double b=pref;
    if (*x>0.0) b+=numer*log(*x);
    if (*x<1.0) b+=(denom-numer*log(1-*x));
    return exp(b);
  }
}

double EfficiencyStatistics::sigma(int i) const {
  if (i==0) return eff0_;
  double seekFor=0;
  double x2=0;
  if (i>0 && i<=6) {
    seekFor=erf(i/sqrt(2));
    x2=1.0;
  } 
  if (i<0 && i>=-6) { 
    seekFor=erf(-i/sqrt(2));
    x2=0.0;
  }
  if (i>0 && numer_==denom_) return 1.0;
  if (i<0 && numer_==0) return 0.0;

  if (denom_>=120) { // use Gaussian-like approximation
    double x;
    if (i<0) x=eff0_+i*0.5*sqrt(eff0_)/sqrt(denom_);
    else x=eff0_+i*0.5*sqrt(1.0-eff0_)/sqrt(denom_);
    if (x<0) x=0;
    if (x>1) x=1;
    return x;
  } else {

    TF1 func("f1",EfficiencyStatFunc,0.0,1.0,2);
    func.SetParameter(0,numer_*1.0);
    func.SetParameter(1,denom_*1.0);
  
    double x1=eff0_;
    const double prec=1e-7;
    double delta,x;
    int cycle=0;
    
    Int_t np = 1000;
    double *lgx=new double[np];
    double *lgw=new double[np];
    func.CalcGaussLegendreSamplingPoints(np,lgx,lgw,1e-15);
    double scale=func.IntegralFast(np,lgx,lgw,x1,x2);
    do {
      x=(x2+x1)/2.0;
      double y=func.IntegralFast(np,lgx,lgw,eff0_,x)/scale;
      delta=fabs(y-seekFor);
      //      printf("%e %f %f %f %f\n",delta,y,x,x1,x2);
      if (y>seekFor) x2=x;
      else x1=x;
      cycle++;
    } while (delta>prec && cycle<20);
    delete [] lgx;
    delete [] lgw;
    return x;
  }
}

double EfficiencyStatistics::rand() const {
  if (denom_==0 || random_==0) return 0;
  double r0,r2;
  double ppeak=prob(eff0_);

  double minv=eff0_-20*rms_;
  if (minv<0) minv=0;
  double maxv=eff0_+20*rms_;
  if (maxv>1) maxv=1;
  double dv=maxv-minv;

  do {    
    r0=random_->Uniform(dv)+minv;
    r2=random_->Uniform(ppeak);
  } while (r2>prob(r0));
  return r0;
}
