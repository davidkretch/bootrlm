#include <algorithm>
#include <cmath>
#include <R.h>
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

/*
 Sampling k from 0:n-1 without replacement.
 */
// copied with modification from MASS/src/lqs.c
// Copyright (C) 1998-2007	B. D. Ripley
// Copyright (C) 1999       R Development Core Team
void sample_noreplace(int *x, int n, int k)
{
  int ind[n];
  int i, j, nn=n;

  for (i = 0; i < n; i++) ind[i] = i;
  for (i = 0; i < k; i++) {
    j = (int)(nn * unif_rand());
    x[i] = ind[j];
    ind[j] = ind[--nn];
  }
}

/*
Find all subsets of size k in order: this gets a new one each call
*/
// copied from MASS/src/lqs.c
// Copyright (C) 1998-2007	B. D. Ripley
// Copyright (C) 1999       R Development Core Team
void next_set(int *x, int n, int k)
{
  int i, j, tmp;

  j = k - 1;
  tmp = x[j]++;
  while(j > 0 && x[j] >= n - (k - 1 -j)) tmp = ++x[--j];
  for(i = j+1; i < k; i++)  x[i] =  ++tmp;
}

} // namespace util
