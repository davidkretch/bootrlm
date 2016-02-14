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
  double qraux[p];
  double work[2*p];
  double res[n];
  double yr[nwhich];
  double xr[nwhich * p];
  double bestcoef[p];
  int pivot[p];
  int which[nwhich];
  //int bestone[nwhich];

  target = (nn - pp)* (beta);

  if(!(sample)) {
    for(i = 0; i < nnew; i++) which[i] = i;
  } else GetRNGstate();

  for(trial = 0; trial < ntrials; trial++) {

    R_CheckUserInterrupt();

    if(!(sample)) {if(trial > 0) util::next_set(which, nn, nnew);}
    else util::sample_noreplace(which, nn, nnew);

    for(j = 0; j < nnew; j++) {
      thisp = which[j];
      yr[j] = y[thisp];
      for(k = 0; k < pp; k++) xr[j + nnew*k] = x[thisp + nn*k];
    }

    /* compute fit, find residuals */
    F77_CALL(dqrdc2)(xr, &nnew, &nnew, &pp, &tol, &rank, qraux, pivot, work);

    if(rank < pp) { sing++; continue; }

    F77_CALL(dqrsl)(xr, &nnew, &nnew, &rank, qraux, yr, &dummy, yr, coef,
             &dummy, &dummy, &n100, &info);

    for(i = 0; i < nn; i++) {
      sum = y[i];
      for(j = 0; j < rank; j++) sum -= coef[j] * x[i + nn*j];
      res[i] = sum;
    }

    /* lqs or lts estimation */
    if(lts < 2) {
      /* find the constant subtracted from the residuals that minimizes
       the criterion. As this is a univariate problem, has an exact
       solution.  */
      if(adj) {
        R_rsort(res, nn);
        if(lts) a = ltsadj(res, nn, qn, &thiscrit);
        else a = lmsadj(res, nn, qn, &thiscrit);
      } else {
        for(i = 0; i < nn; i++) {
          sum = res[i] - a;
          res[i] = sum*sum;
        }
        rPsort(res, nn, qn-1); /* partial sort */
        if(!(lts)) thiscrit = res[qn-1];
        else {
          sum = 0.0;
          for(i = 0; i < qn; i++) sum += res[i];
          thiscrit = sum;
        }
      }
    } else {
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
      } /* now solve for scale S by re-substitution */
      for(iter = 0; iter < 30; iter++) {
        /*printf("iter %d, s = %f sum = %f %f\n", iter, old, sum, target);*/
        sum = 0.0;
        for(i = 0; i < nn; i ++) sum += chi(res[i], k0 * old);
        newp = sqrt(sum/target) * old;
        if(fabs(sum/target - 1.) < 1e-4) break;
        old = newp;
      }
      thiscrit = newp;
    }

    /* first trial might be singular, so use fence */
    if(thiscrit < best) {
      sum = 0.0;
      for(i = 0; i < nn; i ++) sum += chi(res[i], k0 * best);
      best = thiscrit;
      /* printf("trial %d, best = %f sum = %f %f\n", trial, best, sum, target);*/
      //for(i = 0; i < nnew; i++) bestone[i] = which[i] + 1;
      for(i = 0; i < pp; i++) bestcoef[i] = coef[i];
      bestcoef[0] += a;
    }
  } /* for(trial in 0:ntrials) */

  crit = (best < 0.0) ? 0.0 : best;
  if(sample) PutRNGstate();
  /* lqs_free(); */

  // output
  arma::vec coef_out(coef_ptr, p, false);
  arma::vec fitted_out(fitted_ptr, n, false);
  arma::vec resid_out(resid_ptr, n, false);
  arma::vec scale_out(scale_ptr, 1, false);

  for (i = 0; i < p; i++) coef_out[i] = bestcoef[i];
  fitted_out = x * coef_out;
  resid_out = y - fitted_out;
  scale_out = crit;
}

/*
 Adjust the constant for an LMS fit. This is the midpoint of the
 qn contiguous observations of shortest length.
 */
// copied from MASS/src/lqs.c
// Copyright (C) 1998-2007	B. D. Ripley
// Copyright (C) 1999       R Development Core Team
double LQSEstimator::lmsadj(double *x, int n, int qn, double *ssbest)
{
  int i, k = qn - 1;
  double len, best, adj;

  best = x[k] - x[0];
  adj = 0.5*(x[k] + x[0]);
  for(i = 1; i < n-k; i++){
    len = x[i+k] - x[i];
    if(len < best) {
      best = len;
      adj = 0.5*(x[i+k] + x[i]);
    }
  }
  *ssbest = 0.25*best*best;
  return(adj);
}

/*
 Adjust the constant for an LTS fit. This is the mean of the
 qn contiguous observations of smallest variance.
 */
// copied from MASS/src/lqs.c
// Copyright (C) 1998-2007	B. D. Ripley
// Copyright (C) 1999       R Development Core Team
double LQSEstimator::ltsadj(double *x, int n, int qn, double *ssbest)
{
  int i, k = qn - 1;
  double ss, best, m1, m2, adj;

  /*printf("qn = %d\n", qn);*/
  m1 = m2 = 0.0;
  for(i=0; i < qn; i++) {
    m1 += x[i];
    m2 += x[i]*x[i];
  }
  adj = m1/qn;
  best = m2 - m1*m1/qn;

  for(i = 1; i < n-k; i++){
    m1 += x[i+k] - x[i-1];
    m2 += x[i+k]*x[i+k] - x[i-1]*x[i-1];
    ss = m2 - m1*m1/qn;
    if(ss < best) {
      best = ss;
      adj = m1/qn;
    }
  }
  *ssbest = best;
  return(adj);
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
