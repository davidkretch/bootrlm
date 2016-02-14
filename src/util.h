#ifndef UTIL_H
#define UTIL_H

namespace util {

// number of combinations of size k from n elements
unsigned long long choose(unsigned long long, unsigned long long);

// Tukey's biweight function
struct psi {
  psi(double);
  double operator()(double);

private:
  double c;
};

// biweight integral function
struct chi {
  chi(double);
  double operator()(double);

private:
  double c;
};

// sample without replacement
void sample_noreplace(int*, int, int);

// generate all possible combinations
void next_set(int*, int, int);

} // namespace util

#endif
