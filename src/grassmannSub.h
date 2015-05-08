#ifndef GRASSMANSUB_CLASS
#define GRASSMANSUB_CLASS

#include "manifold.h"

//quotient version of grassmann manifold

class grassmannSub: public manifold
{
private:
  arma::mat xiP;
  
public:
  grassmannSub(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the steepest descent direction;
  //double evalObjective(); // evaluate objective
  arma::mat retract(double stepSize,std::string,bool);
  double evalHessian(const arma::mat &H,const arma::mat &Z);
  
  //New added
 // arma::mat genretract(double stepSize, const arma::mat &Z);
  void vectorTrans();
  void update_conjugateD(double);
  void set_particle();
};

#endif
