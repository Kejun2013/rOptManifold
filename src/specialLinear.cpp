#include "specialLinear.h"

specialLinear::specialLinear(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){}


void specialLinear::evalGradient(arma::mat gradF,std::string method){
    double temp=1.0/n*(arma::dot(gradF.t(),Y.i()));
    xi_normal=temp*Y;
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
double specialLinear::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat YU,eucH_proj,Weingarten,Weingarten_part,temp;
  //project euclidian Hessian onto tangent space
  YU=eucH*Y.i();
  YU=1.0/n*(arma::trace(YU))*Y;
  eucH_proj=eucH-YU;
  //Weingarten Map
  Weingarten=-1.0/n*(arma::dot(eucH.t(),Y.i()))*Z;
  Weingarten_part=1.0/n*eucH*Y.i()*Z;
  temp=1.0/n*(arma::dot(Weingarten_part.t(),Y.i()))*Y;
  Weingarten_part-=temp;
  Weingarten+=Weingarten_part;
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}


//second argument unused
arma::mat specialLinear::retract(double stepSize, std::string method, bool first, Function expm){
  if(retraction==1){//Normalization  
    double determ,t;
    t=stepSize; 
    do{
      Yt=Y+t*descD;
      determ=arma::det(Yt);
      t=1.0/2*t;
    }while(determ<0.0);
    Yt=Yt/(pow(determ,1.0/n));
  }else if(retraction==0){//Exponential
    // expm()
  }
  return Yt;
}



//vector transport of conjugate direction, from Y to Yt
//evaluated before accept Y
//project descD onto tangent space at Yt
void specialLinear::vectorTrans(){
  arma::mat temp=descD*Yt.i();
  temp=1.0/n*(arma::trace(temp))*Yt;
  descD=descD-temp;
}


//descD should have been vector transported
//xi is the gradient at new position
void specialLinear::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

double specialLinear::metric(const arma::mat &X1,const arma::mat &X2){
 
 return arma::dot(X1,X2); 
}

void specialLinear::set_particle(){
  arma::mat y_temp=arma::randn(n,p);
  double determ=arma::det(y_temp);
  if(determ<0.0) y_temp.col(1)=-y_temp.col(1);
  Y=y_temp/(pow(determ,1.0/n));
  
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
