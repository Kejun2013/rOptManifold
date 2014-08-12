#ifndef FIXRANK_CLASS
#define FIXRANK_CLASS

#include "manifold.h"

//quotient version of grassmann manifold

class fixRank: public manifold
{
private:
  arma::mat U,V,Ut,Vt;  //Y=U*Sigma*V^T Yt=Ut*Sigma_t*Vt^T
  arma::vec Sigma,Sigma_t;
  arma::mat M,Up,Vp,M_desc,Up_desc,Vp_desc;// decomposition element of tagent vector and descent direction
  arma::mat Qu,Ru,Qv,Rv;// used in retraction
public:
  fixRank(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the steepest descent direction;
  arma::mat retract(double stepSize,std::string,bool first, Function expm);
  double evalHessian(const arma::mat &H,const arma::mat &Z);
  virtual void set_descD(arma::mat );
  void vectorTrans();
  void update_conjugateD(double);
  virtual void retrieve_steepest();
  void set_particle();
  virtual void acceptY();
};

#endif