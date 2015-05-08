#ifndef FIXRANKPSD_CLASS
#define FIXRANKPSD_CLASS

#include "manifold.h"

//quotient version of grassmann manifold

class fixRankPSD: public manifold
{
public:
  fixRankPSD(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the steepest descent direction;
  arma::mat retract(double stepSize,std::string,bool first);
  double evalHessian(const arma::mat &H,const arma::mat &Z);
//  virtual void set_descD(arma::mat );
  void vectorTrans();
  void update_conjugateD(double);
  void set_particle();
//  virtual void acceptY();
};

#endif
