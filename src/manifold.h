#ifndef MANIFOLD_CLASS
#define MANIFOLD_CLASS

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <string>
#include <math>

using namespace Rcpp;
using namespace std;

class manifold
{
protected:
  //Y is the current matrix
  //Yt is the candidate matrix for next step
  arma::mat Y,Yt;  
  //xi and xi_normal: gradient in the tangent space and normal space
  //hessian_Z: hessian operator on Z, as a result of function evalHessian() below
  arma::mat xi,xi_normal,hessian_Z;  
  int n,p,r,retraction;  //Y is a n-by-p matrix of rank r
  double eDescent;//expected descent
public:
  manifold(int n, int p, int r, NumericMatrix,int retraction1);
  
  //evaluate the steepest descent direction;
  virtual void evalGradient(arma::mat gradF, std::string) =0;  
  
  //retraction with certain stepsize
  virtual arma::mat retract(double stepSize) =0;
  
  //evalHessian: Hessian operator on Z at current Y, H is the ordinary hessian
  //            return inner prodcut <Z,Hessian*Z>_Y
  virtual double evalHessian(arma::mat H,arma::mat Z)=0;
  
  //metric on the tagent space at Y: <X1,X2>_Y
  virtual double metric(const arma::mat &X1,const arma::mat &X2)=0; 
  
  //commonly used function across all manifold types
  void acceptY();
  double get_eDescent();
  arma::mat get_hessianZ();
  arma::mat get_Gradient();
  arma::mat get_Y();
  
  //set the descent direction other than the steepest descent direction
  void set_Gradient(arma::mat );
};

#endif