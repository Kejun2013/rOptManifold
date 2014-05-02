#ifndef MANIFOLD_CLASS
#define MANIFOLD_CLASS

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <string>

using namespace Rcpp;
using namespace std;

class manifold
{
protected:
  arma::mat Y,Yt;  //
  arma::mat xi;  //gradient in the tangent space
  int n,p,r,retraction;  //Y is a n-by-p matrix of rank r
  double eDescent;//expected descent
public:
  manifold(int n, int p, int r, NumericMatrix,int retraction1);
  virtual void evalGradient(arma::mat gradF) =0;  //evaluate the steepest descent direction;
  //virtual void evalHessian(arma::mat, arma::mat) =0;//hessian operator
  virtual arma::mat retract(double stepSize) =0;
  void acceptY();
  double get_eDescent();
};

#endif