#include "boot.h"

namespace bootrlm {

// bootstrap class
Boot::Boot(const Data& data, Estimator* estimator, const arma::umat& rep)
  : data(data), n(data.n), p(data.p), rep(rep), r(rep.n_cols) {

  coef.set_size(p, r);
  fitted.set_size(n, r);
  resid.set_size(n, r);
  scale.set_size(1, r);

  arma::uvec indices;

  Data resamp;

  // for each replicate, construct the dataset then estimate
  for (int i = 0; i < r; i++) {
    indices = rep.col(i);
    resamp.x = data.x.rows(indices);
    resamp.y = data.y.elem(indices);
    (*estimator)(resamp, coef.colptr(i), fitted.colptr(i), resid.colptr(i),
                 scale.colptr(i));
  }
}

} // namespace bootrlm