#ifndef DATA_H
#define DATA_H
#include <RcppArmadillo.h>

namespace bootrlm {

class Data {
public:
  arma::vec y;
  arma::mat x;
  int n;
  int p;

  Data();

  Data(arma::vec, arma::mat);
};

} // namespace bootrlm

#endif