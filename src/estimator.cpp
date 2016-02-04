#include "estimator.h"
#include "lqs.h"
#include "util.h"

namespace bootrlm {

// statistical estimator superclass
Estimator::Estimator(std::string type) : type(type) {};

void Estimator::operator()(const Data&, double*, double*, double*, double*) {};

// MM-estimator class
MMEstimator::MMEstimator(const Data& data, double k0)
  : Estimator("MM"), n(data.n), p(data.p), k0(k0) {
  qn = ceil(n/2);
  lts = 2;
  adj = 0;
  sample = (util::choose(n, p) > 5000) ? 1 : 0;
  nwhich = p;
  ntrials = std::min(500 * p, 3000);
  crit = 0;
  sing = 0;
  pk0 = k0;
  beta = 0.5;
}

void MMEstimator::operator()(const Data& data, double* coef_ptr,
                             double* fitted_ptr, double* resid_ptr,
                             double* scale_ptr) {

  const arma::vec& y = data.y;
  const arma::mat& x = data.x;

  arma::vec coef(p);
  arma::vec fitted(n);
  arma::vec resid(n);
  double scale;

  arma::vec y_wt(n);
  arma::mat x_wt(n, p);
  arma::vec wt(n);

  // S estimation
  int bestone[p];
  double bestcoef[p];
  double scale_prev;

  lqs_fitlots(x.memptr(), y.memptr(), &n, &p, &qn, &lts, &adj, &sample, &nwhich,
              &ntrials, &crit, &sing, bestone, bestcoef, &pk0, &beta);

  coef = arma::vec(bestcoef, p);
  resid = y - x * coef;
  scale = crit;

  for (int i = 0; i < 30; i++) {
    // calculate weights
    wt = resid / scale;
    wt.transform(util::psi(k0));
    wt = sqrt(wt);

    // apply weights to data
    y_wt = y % wt;
    // TODO: try sparse diagonal matrix multiplication, column-wise dot product
    for (int i = 0; i < n; i++) {
      x_wt.row(i) = x.row(i) * wt(i);
    }

    // find weighted least squares estimate
    coef = arma::solve(x_wt, y_wt);
    resid = y - x * coef;

    // re-calculate scale
    wt = resid / scale;
    wt.transform(util::chi(k0));
    scale_prev = scale;
    scale *= sqrt(sum(wt) / ((n - p) * beta));

    if (fabs(scale / scale_prev - 1) < 1e-5) break;
  }

  // M estimation
  int maxit = 50000;
  arma::vec resid_prev(n);
  double c = 4.685;
  arma::vec conv(1);

  for (int i = 0; i < maxit; i++) {
    resid_prev = resid;

    // calculate weights
    wt = resid / scale;
    wt.transform(util::psi(c));
    wt = sqrt(wt);

    // apply weights to data
    y_wt = y % wt;
    for (int i = 0; i < n; i++) {
      x_wt.row(i) = x.row(i) * wt(i);
    }

    // find weighted least squares estimate
    coef = arma::solve(x_wt, y_wt);
    resid = y - x * coef;

    // check for convergence
    conv = sqrt(sum(pow(resid - resid_prev, 2)) / sum(pow(resid_prev, 2)));
    if (conv(0) < 1e-4) break;
  }

  // output
  arma::vec coef_out(coef_ptr, p, false);
  arma::vec fitted_out(fitted_ptr, n, false);
  arma::vec resid_out(resid_ptr, n, false);
  arma::vec scale_out(scale_ptr, 1, false);

  coef_out = coef;
  fitted_out = x * coef;
  resid_out = resid;
  scale_out = scale;
}

} // namespace bootrlm