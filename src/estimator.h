#ifndef ESTIMATOR_H
#define ESTIMATOR_H
#include <RcppArmadillo.h>
#include "data.h"

namespace bootrlm {

class Estimator {
public:
  std::string type;
  Estimator(std::string);
  virtual void operator()(const Data&, double*, double*, double*, double*);
};

class MMEstimator : public Estimator {
public:
  MMEstimator(const Data&, double);
  void operator()(const Data&, double*, double*, double*, double*);

private:
  int n;
  int p;
  double k0;
  int qn;
  int lts;
  int adj;
  int sample;
  int nwhich;
  int ntrials;
  double crit;
  int sing;
  double pk0;
  double beta;
};

} // namespace bootrlm

#endif