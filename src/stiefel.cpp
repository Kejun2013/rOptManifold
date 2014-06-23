#include "stiefel.h"

stiefel::stiefel(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){}


void stiefel::evalGradient(arma::mat gradF,std::string method){
    xi_normal=Y.t()*gradF;
    xi_normal=0.5*Y*(xi_normal+xi_normal.t());
    xi=gradF-xi_normal;
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
    }else if(method=="particleSwarm"){
      descD=xi;  //why not -xi?
    }
//    xi=gradF-Y*(gradF.t())*Y;
}


//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double stiefel::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat YU,eucH_proj,Weingarten;
  //project euclidian Hessian onto tangent space
  YU=Y.t()*eucH;
  YU=YU+YU.t();
  eucH_proj=eucH-0.5*Y*YU;
  //Weingarten Map
  YU=Z.t()*xi_normal;
  Weingarten=-Z*Y.t()*xi_normal-0.5*Y*(YU+YU.t());
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}

//second argument unused
arma::mat stiefel::retract(double stepSize, std::string method, bool first){
  if(retraction==1){//QR retraction
    Yt=Y+stepSize*descD;
    arma::qr_econ(retract_Q,retract_R,Yt);
    if(retract_R(0,0)<0){
      retract_Q=-retract_Q;
    }
    Yt=retract_Q;
  }else if(retraction==2){//Cayley
    //Rcpp::Rcout<<"Cayley"<<endl;
    arma::mat Omega=-Y.t()*descD;
    arma::mat temp=arma::eye(p,p)-stepSize/2*Omega;
    Yt=Y*temp.i()*(arma::eye(p,p)+stepSize/2*Omega);
  }
  return Yt;
}

//
//
//
//arma::mat stiefel::genretract(double stepSize, const arma::mat &Z){
//  if(retraction==1){//QR retraction
//    Yt=Y+stepSize*Z;
//    arma::qr_econ(retract_Q,retract_R,Yt);
//    if(retract_R(0,0)<0){
//      retract_Q=-retract_Q;
//    }
//    Yt=retract_Q;
//  }else if(retraction==2){//cayley
//    arma::mat A=Y.t()*Z;  //another approach that z=YA+Y_annhilator*B;
//    arma::mat Omega=(Z-1/2*Y*A)*Y.t()-Y*(Z.t()-1/2*A.t()*Y.t());
//    Yt=arma::eye(n,n)-stepSize/2*Omega;
//    Yt=Yt.i()*(arma::eye(n,n)+stepSize/2*Omega)*Y;
//  }
//  return Yt;
//}


//vector transport of conjugate direction, from Y to Yt
//evaluated before accept Y
//project descD onto tangent space at Yt
void stiefel::vectorTrans(){
  arma::mat temp;
  temp=Yt.t()*descD;
  temp=temp+temp.t();
  descD=descD-0.5*Yt*temp;
}


//descD should have been vector transported
//xi is the gradient at new position
void stiefel::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

double stiefel::metric(const arma::mat &X1,const arma::mat &X2){
// if(retraction==2) return arma::dot(X1,(arma::eye(n,n)-1/2*Y*Y.t())*X2);
// else 
 return arma::dot(X1,X2); 
}

void stiefel::set_particle(){
 arma::mat y_temp=arma::randn(n,p);
  arma::mat Q,R;
  arma::qr_econ(Q,R,y_temp);
  Y=Q;
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
