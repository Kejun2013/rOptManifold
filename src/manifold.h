#ifndef MANIFOLD_CLASS
#define MANIFOLD_CLASS

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <string>
#include <math.h>

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
  //descD: descent direction
  arma::mat xi,xi_normal,hessian_Z,descD;  
  //Y is a n-by-p matrix of rank r
  //retraction is type of retraction expressed in 0,1,2,....
  int n,p,r,retraction;  
  double eDescent,Z_hessian_Z;//expected descent, and <Z, Hessian*Z>_Y
  
public:
  manifold(int n, int p, int r, NumericMatrix,int retraction1);
  
  //evaluate the steepest descent direction;
  virtual void evalGradient(arma::mat gradF, std::string) =0;  
  
  //retraction with certain stepsize
  virtual arma::mat retract(double stepSize, std::string) =0;
  
  //evalHessian: Hessian operator on Z at current Y, H is the ordinary hessian
  //            return inner prodcut <Z,Hessian*Z>_Y
  virtual double evalHessian(const arma::mat &H,const arma::mat &Z)=0;
  
  
  //virtual void vectorTrans()=0;
  
  
  
  //commonly used function across all manifold types
  void acceptY();
  
  //metric on the tagent space at Y: <X1,X2>_Y
  double metric(const arma::mat &X1,const arma::mat &X2); 
  
  
  //second order approximation of the objective function
  double secondOrderApprox(const double m0,const arma::mat &eta,const arma::mat & hessianF);
  
  //return related value or matrix
  double get_eDescent();
  arma::mat get_hessianZ();
  arma::mat get_Gradient();
  arma::mat get_Y();
  
  //set the descent direction other than the steepest descent direction
  void set_descD(arma::mat );
};

#endif