#ifndef ESTIMATOR_H
#define ESTIMATOR_H
#include <RcppArmadillo.h>
#include "data.h"
#include "extern.h"

namespace bootrlm {

class Estimator {
public:
  std::string type;
  Estimator(std::string);
  virtual void operator()(const Data&, double*, double*, double*, double*);
};

class LQSEstimator : public Estimator {
public:
  LQSEstimator();
  LQSEstimator(const Data&, double, double);
  void operator()(const Data&, double*, double*, double*, double*);

private:
  double chi(double, double);
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
  arma::umat indices;
};

class MMEstimator : public Estimator {
public:
  MMEstimator();
  MMEstimator(const Data&, double);
  void operator()(const Data&, double*, double*, double*, double*);

private:
  int n;
  int p;
  double k0;
  double beta;
  LQSEstimator lqs_estimator;
};

} // namespace bootrlm

#endif
