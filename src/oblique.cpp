#include "oblique.h"

oblique::oblique(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){}


void oblique::evalGradient(arma::mat gradF,std::string method){
    arma::mat temp=arma::diagmat(Y.t()*gradF);
    xi_normal=Y*temp;
    xi=gradF-xi_normal;   
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
    }else if(method=="particleSwarm"){
      descD=xi;
    }
//    xi=gradF-Y*(gradF.t())*Y;
}


//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double oblique::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat YU,eucH_proj,Weingarten;
  //project euclidian Hessian onto tangent space
  YU=arma::diagmat(Y.t()*eucH);
  YU=Y*YU;
  eucH_proj=eucH-YU;
  //Weingarten Map
  Weingarten=arma::diagmat(Y.t()*Z);
  Weingarten=-Z*Weingarten;
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}


//second argument unused
arma::mat oblique::retract(double stepSize, std::string method, bool first, Function expm){
  Yt=Y+stepSize*descD;
  arma::mat temp, divisor=arma::eye(p,p);
  int k;
  for(k=0;k<p;k++){
    temp=Yt.col(k);
    divisor(k,k)=arma::norm(temp,"fro");
  }
  Yt=Yt*(divisor.i());
  
  return Yt;
}



//vector transport of conjugate direction, from Y to Yt
//evaluated before accept Y
//project descD onto tangent space at Yt
void oblique::vectorTrans(){
  arma::mat temp=arma::diagmat(Yt.t()*descD);
  temp=Yt*temp;
  descD=descD-temp;
}


//descD should have been vector transported
//xi is the gradient at new position
void oblique::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

double oblique::metric(const arma::mat &X1,const arma::mat &X2){
 
 return arma::dot(X1,X2); 
}

void oblique::set_particle(){
  arma::mat y_temp=arma::randn(n,p), temp, divisor=arma::eye(p,p);
  int k;
  for(k=0;k<p;k++){
    temp=y_temp.col(k);
    divisor(k,k)=arma::norm(temp,"fro");
  }  
  Y=y_temp*(divisor.i());
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
