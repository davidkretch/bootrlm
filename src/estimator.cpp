#include "estimator.h"
#include "util.h"

namespace bootrlm {

// statistical estimator superclass
Estimator::Estimator(std::string type) : type(type) {};

void Estimator::operator()(const Data&, double*, double*, double*, double*) {};

// LQS-estimator class
LQSEstimator::LQSEstimator() : Estimator("LQS") {};

LQSEstimator::LQSEstimator(const Data& data, double k0, double beta)
  : Estimator("LQS"), n(data.n), p(data.p), k0(k0), beta(beta) {
  qn = ceil(n/2);
  lts = 2;
  adj = 0;
  sample = (util::choose(n, p) > 5000) ? 1 : 0;
  nwhich = p;
  ntrials = std::min(500 * p, 3000);
  crit = 0;
  sing = 0;
  pk0 = k0;

  indices.set_size(nwhich, ntrials);
  int which[nwhich];
  arma::uvec which_vec(nwhich);

  if(!(sample)) {
    for(int i = 0; i < nwhich; i++) which[i] = i;
  } else GetRNGstate();

  // get all the subset indices used across all trials
  for (int trial = 0; trial < ntrials; trial++) {
    if(!(sample)) {if(trial > 0) util::next_set(which, n, nwhich);}
    else util::sample_noreplace(which, n, nwhich);

    for (int i = 0; i < nwhich; i++) which_vec[i] = (unsigned int) which[i];
    indices.col(trial) = which_vec;
  }
}

// For lots of subsets of size nwhich, compute the exact fit to those data
// points and the residuals from all the data points.
// copied with modification from MASS/src/lqs.c
// Copyright (C) 1998-2007	B. D. Ripley
// Copyright (C) 1999       R Development Core Team
// TODO: rewrite
void LQSEstimator::operator()(const Data& data, double* coef_ptr,
                            double* fitted_ptr, double* resid_ptr,
                            double* scale_ptr) {
  int nnew = nwhich, pp = p;
  int i, iter, j, k,  nn = n, thisp, trial;
  int rank, info, n100 = 100;
  int firsttrial = 1;
  double a = 0.0, tol = 1.0e-7, sum, thiscrit, best = DBL_MAX, target, old,
    newp, dummy, k0 = pk0;

  const arma::vec& y = data.y;
  const arma::mat& x = data.x;

  double coef[p];
  arma::vec coef_vec(coef, p, false, true);
  double qraux[p];
  double work[2*p];
  double res[n];
  arma::vec res_vec(res, n, false, true);
  double yr[nwhich];
  double xr[nwhich * p];
  arma::vec yr_vec(yr, nwhich, false, true);
  arma::mat xr_mat(xr, nwhich, p, false, true);
  double bestcoef[p];
  int pivot[p];
  arma::uvec which_vec(nwhich);
  //int bestone[nwhich];

  target = (nn - pp)* (beta);

  for(trial = 0; trial < ntrials; trial++) {

    R_CheckUserInterrupt();

    // get this trial's subset
    which_vec = indices.col(trial);
    yr_vec = y.elem(which_vec);
    xr_mat = x.rows(which_vec);

    /* compute fit, find residuals */
    F77_CALL(dqrdc2)(xr, &nnew, &nnew, &pp, &tol, &rank, qraux, pivot, work);

    if(rank < pp) { sing++; continue; }

    F77_CALL(dqrsl)(xr, &nnew, &nnew, &rank, qraux, yr, &dummy, yr, coef,
             &dummy, &dummy, &n100, &info);

    res_vec = y - x * coef_vec;

    /* S estimation */
    if(firsttrial) {
      for(i = 0; i < nn; i ++) res[i] = fabs(res[i]);
      rPsort(res, nn, nn/2);
      old = res[nn/2]/0.6745;	 /* MAD provides the initial scale */
      firsttrial = 0;
    } else {
      /* only find optimal scale if it will be better than
       existing best solution */
      sum = 0.0;
      for(i = 0; i < nn; i ++) sum += chi(res[i], k0 * best);
      if(sum > target) continue;
      old = best;
    }

    /* now solve for scale S by re-substitution */
    for(iter = 0; iter < 30; iter++) {
      /*printf("iter %d, s = %f sum = %f %f\n", iter, old, sum, target);*/
      sum = 0.0;
      for(i = 0; i < nn; i ++) sum += chi(res[i], k0 * old);
      newp = sqrt(sum/target) * old;
      if(fabs(sum/target - 1.) < 1e-4) break;
      old = newp;
    }
    thiscrit = newp;

    /* first trial might be singular, so use fence */
    if(thiscrit < best) {
      sum = 0.0;
      for(i = 0; i < nn; i ++) sum += chi(res[i], k0 * best);
      best = thiscrit;
      for(i = 0; i < pp; i++) bestcoef[i] = coef[i];
      bestcoef[0] += a;
    }
  } /* for(trial in 0:ntrials) */

  crit = (best < 0.0) ? 0.0 : best;
  if(sample) PutRNGstate();
  /* lqs_free(); */

  // output
  arma::vec coef_out(coef_ptr, p, false, true);
  arma::vec fitted_out(fitted_ptr, n, false, true);
  arma::vec resid_out(resid_ptr, n, false, true);
  arma::vec scale_out(scale_ptr, 1, false, true);

  for (i = 0; i < p; i++) coef_out[i] = bestcoef[i];
  fitted_out = x * coef_out;
  resid_out = y - fitted_out;
  scale_out = crit;
}

/* the chi function for the S estimator: the integral of biweight */
// copied from MASS/src/lqs.c
// Copyright (C) 1998-2007	B. D. Ripley
// Copyright (C) 1999       R Development Core Team
// TODO: find a more appropriate location
double LQSEstimator::chi(double x, double a)
{
  x /= a; x *= x;
  if(x > 1) return(1.0);
  else return(x*(3 + x*(-3 + x)));
}

// MM-estimator class
MMEstimator::MMEstimator() : Estimator("MM") {}

MMEstimator::MMEstimator(const Data& data, double k0)
  : Estimator("MM"), n(data.n), p(data.p), k0(k0) {
  beta = 0.5;
  lqs_estimator = LQSEstimator(data, k0, beta);
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
  double scale_prev;

  lqs_estimator(data, coef.memptr(), fitted.memptr(), resid.memptr(), &scale);

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
  arma::vec coef_out(coef_ptr, p, false, true);
  arma::vec fitted_out(fitted_ptr, n, false, true);
  arma::vec resid_out(resid_ptr, n, false, true);
  arma::vec scale_out(scale_ptr, 1, false, true);

  coef_out = coef;
  fitted_out = x * coef;
  resid_out = resid;
  scale_out = scale;
}

} // namespace bootrlm
