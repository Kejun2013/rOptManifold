#ifndef FIXRANKSYM_CLASS
#define FIXRANKSYM_CLASS

#include "manifold.h"

//quotient version of grassmann manifold

class fixRankSym: public manifold
{
private:
  arma::mat V,Vt;  //Y=U*Sigma*V^T Yt=Ut*Sigma_t*Vt^T
  arma::vec Sigma,Sigma_t;
  arma::mat M,Vp,M_desc,Vp_desc;// decomposition element of tagent vector and descent direction
  arma::mat Qv,Rv;// used in retraction
public:
  fixRankSym(int, int, int, NumericMatrix,int);
  void evalGradient(arma::mat gradF,std::string method);  //evaluate the steepest descent direction;
  arma::mat retract(double stepSize,std::string,bool first);
  double evalHessian(const arma::mat &H,const arma::mat &Z);
  virtual void set_descD(arma::mat );
  void vectorTrans();
  void update_conjugateD(double);
  virtual void retrieve_steepest();
  void set_particle();
  virtual void acceptY();
};

#endif