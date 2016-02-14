#ifndef EXTERN_H
#define EXTERN_H

namespace bootrlm {

extern "C" {

  // QR matrix decomposition functions
  void F77_NAME(dqrsl)(double*, int*, int*, int*, double*, double*, double*,
                double*, double*, double*, double*, int*, int*);

  void F77_NAME(dqrdc2)(double*, int*, int*, int*, double*, int*, double*, int*,
                double*);

} // extern "C"

} // namespace bootrlm

#endif
