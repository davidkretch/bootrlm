#include "boot.h"

// TODO: clean up interface to output matrices

namespace bootrlm {

// bootstrap class
Boot::Boot(const Data& data, Estimator* estimator, const arma::umat& rep)
  : data(data), n(data.n), p(data.p), estimator(estimator), rep(rep), r(rep.n_cols) {

  coef.set_size(p, r);
  fitted.set_size(n, r);
  resid.set_size(n, r);
  scale.set_size(1, r);

  ProcessReplicate process_replicate(data, estimator, rep, coef, fitted, resid, scale);

  // for each replicate, construct the dataset then estimate
  RcppParallel::parallelFor(0, r, process_replicate);

}

Boot::ProcessReplicate::ProcessReplicate(const Data& data,
                                         Estimator* estimator,
                                         const arma::umat& rep,
                                         arma::mat& coef,
                                         arma::mat& fitted,
                                         arma::mat& resid,
                                         arma::mat& scale)
  : data(data), estimator(estimator), rep(rep), coef(coef), fitted(fitted),
    resid(resid), scale(scale) {};

void Boot::ProcessReplicate::operator()(std::size_t begin, std::size_t end) {
  arma::uvec indices;
  Data replicate;

  for (std::size_t i = begin; i < end; ++i) {
    indices = rep.col(i);
    replicate.x = data.x.rows(indices);
    replicate.y = data.y.elem(indices);
    (*estimator)(replicate, coef.colptr(i), fitted.colptr(i),
     resid.colptr(i), scale.colptr(i));
  }
}

} // namespace bootrlm
