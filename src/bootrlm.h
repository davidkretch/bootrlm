#ifndef BOOTRLM_H
#define BOOTRLM_H
#include <RcppArmadillo.h>


class Data {
public:
  arma::vec y;
  arma::mat x;
  int n;
  int p;

  Data();

  Data(arma::vec, arma::mat);
};


class Estimator {
public:
  std::string type;
  Estimator(std::string);
  virtual void operator()(const Data&, double*, double*, double*, double*);
};


class MMEstimator : public Estimator {
public:
  MMEstimator(const Data&, double);
  void operator()(const Data&, double*, double*, double*, double*);

private:
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
};


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


Rcpp::List bootrlm_cpp(arma::vec, arma::mat, arma::umat, int, Rcpp::List);


#endif
