#ifndef BOOTRLM_H
#define BOOTRLM_H
#include <RcppArmadillo.h>

Rcpp::List bootrlm_cpp(arma::vec, arma::mat, arma::umat, int, Rcpp::List);

#endif
