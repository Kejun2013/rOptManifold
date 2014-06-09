#include "grassmannQ.h"

void grassmannQ::evalGradient(arma::mat gradF, std::string method){
  arma::mat U_svd;
  xi=gradF-Y*Y.t()*gradF;
  eDescent=arma::dot(gradF,xi);
  arma::svd_econ(U_svd,retract_Sigma,retract_V,-xi);
  retract_U=arma::mat(n,2*p,arma::fill::zeros);
  retract_U(arma::span::all,arma::span(0,p-1))=Y*retract_V;
  retract_U(arma::span::all,arma::span(p,2*p-1))=U_svd;
}



arma::mat grassmannQ::retract(double stepSize, std::string method,bool first){
  arma::mat retract_middle;
  retract_middle=arma::mat(2*p,p,arma::fill::zeros);
  retract_middle(arma::span(0,p-1),arma::span::all)
    =diagmat(cos(retract_Sigma*stepSize));
  retract_middle(arma::span(p,2*p-1),arma::span::all)
    =diagmat(sin(retract_Sigma*stepSize));
  Yt=retract_U*retract_middle*retract_V.t();
  return Yt;
}

grassmannQ::grassmannQ(int n1, int p1, int r1, 
          Rcpp::NumericMatrix Y1,int retraction1):
          manifold(n1,p1,r1,Y1,retraction1){}


double grassmannQ::evalHessian(const arma::mat &H, const arma::mat &Z){
  // something here
  return 0;
}


arma::mat grassmannQ::genretract(double stepSize, const arma::mat &Z){
//empty
  return arma::eye(n,n);
}

void grassmannQ::vectorTrans(){
    //empty
}

void grassmannQ::update_conjugateD(double eta){
  //empty
}

void grassmannQ::set_particle(){
  Yt=Y;
}

