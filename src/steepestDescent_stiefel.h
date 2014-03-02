#ifndef _rOptManifold_Stiefel_H
#define _rOptManifold_Stiefel_H

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>

RcppExport SEXP steepestDescent_stiefel(SEXP n1, SEXP p1, SEXP f1, SEXP f2, SEXP f3,
                          SEXP iterMax1, SEXP tol1, SEXP retract1);

#endif
