#ifndef LQS_H
#define LQS_H

void
  lqs_fitlots(const double *x, const double *y, int *n, int *p, int *qn,
              int *lts, int *adj, int *sample, int *nwhich,
              int *ntrials, double *crit, int *sing, int *bestone,
              double *bestcoef, double *pk0, double *beta);

double
  chi(double x, double a);

#endif
