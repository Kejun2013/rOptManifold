#include "manifold.h"

//m_x(eta)=m0+<xi,eta>+0.5*<eta,Hessian*eta>
double manifold::secondOrderApprox(const double m0,const arma::mat &eta,
                                    const arma::mat & hessianF){
  evalHessian(hessianF,eta);
  double f;
  f=m0+0.5*Z_hessian_Z+metric(eta,xi);
  return f;
}

double manifold::metric(const arma::mat &X1,const arma::mat &X2){
  return arma::dot(X1,X2);
}


manifold::manifold(int n1, int p1, int r1, 
              Rcpp::NumericMatrix Y1,int retraction1){
  n=n1;
  p=p1;
  r=r1;
  Y=arma::mat(Y1.begin(),n,p,false);
  retraction=retraction1;
}


void manifold::acceptY(){
  Y=Yt;
}


arma::mat manifold::get_Y(){
  return Y;
}


double manifold::get_eDescent(){
  return eDescent;
}


arma::mat manifold::get_Gradient(){
  return xi;
}

arma::mat manifold::get_hessianZ(){
  return hessian_Z;
}


//set descent direction
void manifold::set_descD(arma::mat xi_temp){
  descD=xi_temp;
}


