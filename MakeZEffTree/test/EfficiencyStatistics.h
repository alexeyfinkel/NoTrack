/* -*- C++ -*- */
#ifndef EfficiencyStatistics_h_included
#define EfficiencyStatistics_h_included 1

#include <vector>
class TRandom;

class EfficiencyStatistics {
public:
  EfficiencyStatistics(double eff0, int denom);
  EfficiencyStatistics(int numer, int denom);
  
  void setRandom(TRandom* rand);
  double mean() const { return eff0_; }
  double sigma(int i) const;
  double rand() const;
  double prob(double eff) const;
private:
  double eff0_,rms_;
  double prefactor_;
  int numer_, denom_;
  TRandom* random_;
};

#endif // EfficiencyStatistics_h_included
