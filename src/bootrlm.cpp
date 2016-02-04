//' @useDynLib bootrlm
//' @importFrom Rcpp sourceCpp
#include "boot.h"
#include "data.h"
#include "estimator.h"
#include "bootrlm.h"

using namespace bootrlm;

// main function and interface to R
// [[Rcpp::export]]
Rcpp::List bootrlm_cpp(arma::vec y, arma::mat x, arma::umat rep, int method,
                       Rcpp::List options) {

  Data data = Data(y, x);

  Estimator* estimator = NULL;

  switch(method) {
  case 1:
    estimator = new MMEstimator(data, options["k"]);
  }

  Boot boot = Boot(data, estimator, rep);

  delete estimator;

  Rcpp::List ret;
  ret["r"] = boot.r;
  ret["replicates"] = boot.rep;
  ret["coefficients"] = boot.coef;
  ret["fitted"] = boot.fitted;
  ret["residuals"] = boot.resid;
  ret["scale"] = boot.scale;

  return ret;
};