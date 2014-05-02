#ifndef GRASSMANQ_CLASS
#define GRASSMANQ_CLASS

#include <manifold.h>

class grassmanQ: public manifold
{
private:
  arma::mat retract_U,retract_V;
  arma::vec retract_Sigma;  //gradient in the tangent space
public:
  grassmanQ(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF);  //evaluate the steepest descent direction;
  //double evalObjective(); // evaluate objective
  arma::mat retract(double stepSize);
};

#endif