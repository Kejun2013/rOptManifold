#include "sphere.h"

sphere::sphere(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){}


void sphere::evalGradient(arma::mat gradF,std::string method){
    xi_normal=(arma::dot(Y,gradF))*Y;
    xi=gradF-xi_normal;
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
    }else if(method=="particleSwarm"){
      descD=-xi;
    }
//    xi=gradF-Y*(gradF.t())*Y;
}


//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double sphere::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat YU,eucH_proj,Weingarten;
  //project euclidian Hessian onto tangent space
  YU=(arma::dot(Y,eucH))*Y;
  eucH_proj=eucH-YU;
  //Weingarten Map
  Weingarten=-(arma::dot(Y,Z))*Z;
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}


//second argument unused
arma::mat sphere::retract(double stepSize, std::string method, bool first){
  if(retraction==1){//Normalization  //?How to change in setMethod, add "nor=1"?
    Yt=Y+stepSize*descD;
    double norm;
    norm=arma::dot(Yt,Yt);
    norm=pow(norm,1/2.0);
    Yt=Yt/norm;
  }else if(retraction==0){//Exponential
    double norm;
    norm=arma::dot(stepSize*descD,stepSize*descD);
    norm=pow(norm,1/2.0);
    double norm_1=arma::dot(descD,descD);
    norm_1=pow(norm_1,1/2.0);    
    Yt=cos(norm)*Y+sin(norm)*descD/norm_1;
  }
  return Yt;
}




arma::mat sphere::genretract(double stepSize, const arma::mat &Z){
  if(retraction==1){//QR retraction
    Yt=Y+stepSize*Z;
    double norm;
    norm=arma::dot(Yt,Yt);
    norm=pow(norm,1/2.0);
    Yt=Yt/norm;
  }else if(retraction==0){//Exponential
    double norm;
    norm=arma::dot(stepSize*Z,stepSize*Z);
    norm=pow(norm,1/2.0);
    double norm_1=arma::dot(Z,Z);
    norm_1=pow(norm_1,1/2.0);
    Yt=cos(norm)*Y+sin(norm)*Z/norm_1;
  }
  return Yt;
}


//vector transport of conjugate direction, from Y to Yt
//evaluated before accept Y
//project descD onto tangent space at Yt
void sphere::vectorTrans(){
  arma::mat temp;
  temp=(arma::dot(Yt,descD))*Yt;
  descD=descD-temp;
}


//descD should have been vector transported
//xi is the gradient at new position
void sphere::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

double sphere::metric(const arma::mat &X1,const arma::mat &X2){
// if(retraction==2) return arma::dot(X1,(arma::eye(n,n)-1/2*Y*Y.t())*X2);
// else 
 return arma::dot(X1,X2); 
}

void sphere::set_particle(){
  arma::mat y_temp=arma::randn(n,p);
  double norm;
  norm=arma::dot(y_temp,y_temp);
  norm=pow(norm,1/2.0);
  Y=y_temp/norm;
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
