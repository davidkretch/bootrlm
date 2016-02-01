#include <algorithm>
#include <cmath>
#include "util.h"


namespace util {


// number of combinations of size k from n elements
unsigned long long choose(unsigned long long n, unsigned long long k) {
  if (k > n) {
    return 0;
  }
  unsigned long long r = 1;
  for (unsigned long long d = 1; d <= k; ++d) {
    r *= n--;
    r /= d;
  }
  return r;
}


// initialize c for Tukey's biweight function
psi::psi(double c) : c(c) {};


// Tukey's biweight function
double psi::operator()(double x) {
  return pow(1 - pow(std::min(1.0, fabs(x / c)), 2), 2);
};


// initialize c for biweight integral function
chi::chi(double c) : c(c) {};


// biweight integral function
double chi::operator()(double x) {
  x /= c;
  x *= x;
  if (x > 1) return 1.0;
  else return x * (3 + x * (-3 + x));
};


}
