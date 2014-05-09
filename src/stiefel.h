#ifndef STIEFEL_CLASS
#define STIEFEL_CLASS

#include <manifold.h>

class stiefel: public manifold
{
private:
  arma::mat retract_Q,retract_R;  //gradient in the tangent space
public:
  stiefel(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF);  //evaluate the steepest descent direction;
  //double evalObjective(); // evaluate objective
  arma::mat retract(double stepSize);
  double evalHessian(arma::mat H,arma::mat Z);
  double metric(const arma::mat &,const arma::mat &);
};

#endif