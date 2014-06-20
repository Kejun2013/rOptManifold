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
  //descD: descent direction, used as conjugate direction in conjugate descent algorithm
  //conjugdateD: conjugate descent direction
  arma::mat xi,xi_normal,hessian_Z,descD;
  //Y is a n-by-p matrix of rank r
  //retraction is type of retraction expressed in 0,1,2,3,....
  //retraction: 
  int n,p,r,retraction;  
  double eDescent,Z_hessian_Z;//expected descent, and <Z, Hessian*Z>_Y
  
public:
  manifold(int n, int p, int r, NumericMatrix ,int retraction1);
  
  //evaluate the steepest descent direction;
  virtual void evalGradient(arma::mat gradF, std::string) =0;  
  
  //retraction with certain stepsize
  virtual arma::mat retract(double stepSize, std::string,bool) =0;
  
  //evalHessian: Hessian operator on Z at current Y, H is the ordinary hessian
  //            return inner prodcut <Z,Hessian*Z>_Y
  virtual double evalHessian(const arma::mat &H,const arma::mat &Z)=0;

 //'general retraction with certain stepsize
 //'
 //'@param stepSize
 //'@param direction of vector
 // virtual arma::mat genretract(double stepSize, const arma::mat &Z)=0;

//' The vector Transportation from Current Y to Yt
//'
//' Need to use genretract() member function
//' @param stepSize
//' @param direction of vector trasportation 
//' @return a tangent vector to the Yt associated with retraction 
  virtual void vectorTrans()=0;
  
//' Set the conjugate direction
//' 
//' @param set eta as the new conjugate (tangent) direction conjugate_temp
  virtual void update_conjugateD (double eta)=0;
  virtual void retrieve_steepest(){descD=-xi;}
//' Set the eDescent
//' 
//' @param set eDescent_temp as the new eDescent
  void set_eDescent (double eDescent_temp) {eDescent=eDescent_temp;}

//set the descent direction other than the steepest descent direction
 virtual  void set_descD(arma::mat );
    
//' Set initial particle swarm;
  virtual void set_particle ()=0;    
  
  void set_Y(arma::mat Y_temp){Y=Y_temp;}
  

virtual  void acceptY();
  
  //metric on the tagent space at Y: <X1,X2>_Y
  //'Kejun puts it to virtual
virtual  double metric(const arma::mat &X1,const arma::mat &X2); 
  double gradMetric(){return metric(xi,xi);}
  double descDMetric(){return metric(descD,descD);}
  double grad_descD_Metric(){return metric(xi,descD);}
  //second order approximation of the objective function
  double secondOrderApprox(const double m0,const arma::mat &eta,const arma::mat & hessianF);
  
  //return related value or matrix
  double get_eDescent();
  arma::mat get_hessianZ();
  arma::mat get_Gradient();
  arma::mat get_Y();
  arma::mat get_descD(){return descD;}
  

};

#endif