#ifndef BOOT_H
#define BOOT_H
#include <RcppArmadillo.h>
#include "data.h"
#include "estimator.h"

namespace bootrlm {

class Boot {
public:
  const Data data;
  int n;
  int p;
  const arma::umat rep;
  int r;
  arma::mat coef;
  arma::mat fitted;
  arma::mat resid;
  arma::mat scale;

  Boot(const Data&, Estimator*, const arma::umat&);
};

} // namespace bootrlm

#endif