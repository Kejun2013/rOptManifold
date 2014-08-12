#ifndef SPHERE_CLASS
#define SPHERE_CLASS

#include "manifold.h"

class sphere: public manifold
{
private:
  // arma::mat retract_Q,retract_R;  //gradient in the tangent space
public:
  sphere(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the steepest descent direction;
  //double evalObjective(); // evaluate objective
  arma::mat retract(double stepSize,std::string method, bool first, Function expm);
  double evalHessian(const arma::mat &H,const arma::mat &Z);

  // new added
  arma::mat genretract(double stepSize, const arma::mat &Z);
  void vectorTrans();
  void update_conjugateD(double);
  double metric(const arma::mat &X1,const arma::mat &X2);
  void set_particle();
};

#endif