#ifndef SPECIALLINEAR_CLASS
#define SPECIALLINEAR_CLASS

#include <manifold.h>

class specialLinear: public manifold
{

public:
  //n*n matrix  p=n has been set up "R" files
  specialLinear(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the tangent vecotors;
  //double evalObjective(); // evaluate objective
  arma::mat retract(double stepSize,std::string method, bool first);
  double evalHessian(const arma::mat &H,const arma::mat &Z);

  // new added
  arma::mat genretract(double stepSize, const arma::mat &Z);
  void vectorTrans();
  void update_conjugateD(double);
  double metric(const arma::mat &X1,const arma::mat &X2);
  void set_particle();
};

#endif

