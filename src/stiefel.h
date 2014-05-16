#ifndef STIEFEL_CLASS
#define STIEFEL_CLASS

#include <manifold.h>

class stiefel: public manifold
{
private:
  arma::mat retract_Q,retract_R;  //gradient in the tangent space
public:
  stiefel(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the steepest descent direction;
  //double evalObjective(); // evaluate objective
  arma::mat retract(double stepSize,std::string method);
  double evalHessian(const arma::mat &H,const arma::mat &Z);
  
  
  // new added
  arma::mat genretract(double stepSize, const arma::mat &Z);
  arma::mat vectorTrans(double stepSize, const arma::mat &Z);
};

#endif