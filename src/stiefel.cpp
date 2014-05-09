#include "stiefel.h"

void stiefel::evalGradient(arma::mat gradF,std::string method){
  
  //xi=xi;
  if(method=="steepest"){
    xi=gradF-Y*gradF.t()*Y;
    eDescent=arma::dot(gradF,xi);
    if(retraction==2){//prepare for Cayley rectraction in steepest descent algorithm
      retract_Q=arma::mat(n,2*p,arma::fill::zeros);
      retract_R=retract_Q;
      retract_Q(arma::span::all,arma::span(0,p-1))=gradF;
      retract_Q(arma::span::all,arma::span(p,2*p-1))=Y;
      retract_R(arma::span::all,arma::span(0,p-1))=Y;
      retract_R(arma::span::all,arma::span(p,2*p-1))=-gradF;
    }
  }else if(method=="trustRegion"){
    //xi_normal=0.5*Y(gradF^T*Y+Y^T*gradF)
    xi_normal=Y.t()*gradF;
    xi_normal=0.5*Y*(xi_normal+xi_normal.t());
    xi=gradF-xi_normal;
  }
}


//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double stiefel::evalHessian(arma::mat eucH,arma::mat Z){
  arma::mat YU,eucH_proj,Weingarten;
  //project euclidian Hessian onto tangent space
  YU=Y.t()*eucH;
  YU=YU+YU.t();
  eucH_proj=eucH-0.5*Y*YU;
  //Weigarten Map
  YU=Z.t()*xi_normal;
  Weingarten=-Z*Y.t()*xi_normal-0.5*Y*(YU+YU.t());
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
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
{}


double arma::metric(const arma::mat &X1, const arma::mat &X2){
  return arma::dot(X1,X2)
}

