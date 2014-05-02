#include "stiefel.h"

void stiefel::evalGradient(arma::mat gradF){
  xi=gradF-Y*gradF.t()*Y;
  //xi=xi;
  eDescent=arma::dot(gradF,xi);
  if(retraction==2){
    retract_Q=arma::mat(n,2*p,arma::fill::zeros);
    retract_R=retract_Q;
    retract_Q(arma::span::all,arma::span(0,p-1))=gradF;
    retract_Q(arma::span::all,arma::span(p,2*p-1))=Y;
    retract_R(arma::span::all,arma::span(0,p-1))=Y;
    retract_R(arma::span::all,arma::span(p,2*p-1))=-gradF;
  }
}



arma::mat stiefel::retract(double stepSize){
  if(retraction==1){//QR retraction
    Yt=Y-stepSize*xi;
    arma::qr_econ(retract_Q,retract_R,Yt);
    if(retract_R(0,0)<0){
      retract_Q=-retract_Q;
    }
    Yt=retract_Q;
  }else if(retraction==2){//Cayley Retraction
    Yt=arma::eye<arma::mat>(2*p,2*p)+stepSize/2*retract_R.t()*retract_Q;
    Yt=Y-stepSize*retract_Q*Yt.i()*retract_R.t()*Y;
  }
  return Yt;
}

stiefel::stiefel(int n1, int p1, int r1, 
          Rcpp::NumericMatrix Y1,int retraction1):
          manifold(n1,p1,r1,Y1,retraction1)
{    
  
}


