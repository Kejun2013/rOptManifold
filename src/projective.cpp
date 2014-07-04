#include "projective.h"

projective::projective(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){}


void projective::evalGradient(arma::mat gradF,std::string method){
    arma::mat temp=(Y.t()*Y);
    temp=temp.i();
    temp=temp.t();
    xi=1.0/2*Y*temp*(gradF.t()*Y+Y.t()*gradF);
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
double projective::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  return 1.0;
  
}


//second argument unused
arma::mat projective::retract(double stepSize, std::string method, bool first){
  Yt=Y+stepSize*descD;

  return Yt;
}



//vector transport of conjugate direction, from Y to Yt
//evaluated before accept Y
//project descD onto tangent space at Yt
void projective::vectorTrans(){
  arma::mat temp=(Yt.t()*Yt);
  temp=temp.i();
  temp=temp.t();
  descD=1.0/2*Yt*temp*(descD.t()*Yt+Yt.t()*descD);
}


//descD should have been vector transported
//xi is the gradient at new position
void projective::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

double projective::metric(const arma::mat &X1,const arma::mat &X2){
 
 return arma::dot(X1,X2); 
}

void projective::set_particle(){
  arma::mat y_temp=arma::randn(n,p);
  
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
