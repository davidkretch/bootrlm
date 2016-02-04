#include "data.h"

namespace bootrlm {

// data container class
Data::Data() {};

Data::Data(arma::vec y, arma::mat x) : y(y), x(x) {
  n = x.n_rows;
  p = x.n_cols;
};

} // namespace bootrlm