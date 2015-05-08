#include "projective.h"

projective::projective(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){}

arma::mat projective::get_horizontal(const arma::mat &H){
  arma::mat H_horizontal,temp=Y.t()*Y;
  H_horizontal=temp.i()*Y.t()*H+H.t()*Y*temp.i();
  H_horizontal=1/2.0*Y*H_horizontal;
  H_horizontal+=(arma::eye(n,n)-Y*temp.i()*Y.t())*H;
  return H_horizontal;
}

void projective::evalGradient(arma::mat gradF,std::string method){
    xi=get_horizontal(gradF);
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
    }else if(method=="particleSwarm"){
      descD=xi;
    }
}


//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double projective::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat Z1,temp=Y.t()*Y;
  temp=(1/2.0)*Z*(temp.i()*Y.t()*eucH-eucH.t()*Y*temp.i());
  temp=-eucH+temp;
  hessian_Z=get_horizontal(temp);
  Z_hessian_Z=metric(Z,hessian_Z);
  return Z_hessian_Z;
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
  arma::mat temp=Yt.t()*Yt,H_horizontal;
  H_horizontal=temp.i()*Yt.t()*descD+descD.t()*Yt*temp.i();
  H_horizontal=1/2.0*Yt*H_horizontal;
  H_horizontal+=(arma::eye(n,n)-Yt*temp.i()*Yt.t())*descD;
  descD=H_horizontal;
}


//descD should have been vector transported
//xi is the gradient at new position
void projective::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

double projective::metric(const arma::mat &X1,const arma::mat &X2){
 arma::mat temp1,temp2,temp=Y.t()*Y;
 temp1=temp.i()*Y.t()*X1;
 temp2=temp.i()*Y.t()*X2;
 return arma::dot(temp1,temp2);
}

void projective::set_particle(){
  arma::mat y_temp=arma::randn(n,p);
  
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
