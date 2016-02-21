#ifndef BOOT_H
#define BOOT_H
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "data.h"
#include "estimator.h"

namespace bootrlm {

class Boot {
public:
  const Data data;
  int n;
  int p;
  Estimator* estimator;
  const arma::umat rep;
  int r;
  arma::mat coef;
  arma::mat fitted;
  arma::mat resid;
  arma::mat scale;

  Boot(const Data&, Estimator*, const arma::umat&);

private:
  struct ProcessReplicate : public RcppParallel::Worker {

    // need to pass class members to struct, as they're not visible (in C++98)
    ProcessReplicate(const Data&, Estimator*, const arma::umat&,
                     arma::mat&, arma::mat&, arma::mat&, arma::mat&);

    void operator()(std::size_t begin, std::size_t end);

    const Data& data;
    Estimator* estimator;
    const arma::umat& rep;
    arma::mat& coef;
    arma::mat& fitted;
    arma::mat& resid;
    arma::mat& scale;
  };
};

} // namespace bootrlm

#endif
