#include "grassmannQ.h"

void grassmannQ::evalGradient(arma::mat gradF, std::string method){
  xi=gradF-Y*Y.t()*gradF;
  eDescent=arma::dot(gradF,xi);
  if(method=="steepest"){
    descD=-xi;
  }else if(method=="particleSwarm"){
    descD=xi;
  }else if(method=="trustRegion"){
    
  }
}



arma::mat grassmannQ::retract(double stepSize, std::string method,bool first){
  if(retraction==1){ //QR
        Yt=Y+stepSize*descD;
    arma::qr_econ(retract_U,retract_V,Yt);
    if(retract_V(0,0)<0){
      retract_U=-retract_U;
    }
    Yt=retract_U;
  }else if(retraction==0){
    if(first){
      arma::mat U_svd;
      arma::svd_econ(U_svd,retract_Sigma,retract_V,-xi);
      retract_U=arma::mat(n,2*p,arma::fill::zeros);
      retract_U(arma::span::all,arma::span(0,p-1))=Y*retract_V;
      retract_U(arma::span::all,arma::span(p,2*p-1))=U_svd;
    }
      arma::mat retract_middle;
      retract_middle=arma::mat(2*p,p,arma::fill::zeros);
      retract_middle(arma::span(0,p-1),arma::span::all)
        =diagmat(cos(retract_Sigma*stepSize));
      retract_middle(arma::span(p,2*p-1),arma::span::all)
        =diagmat(sin(retract_Sigma*stepSize));
      Yt=retract_U*retract_middle*retract_V.t();
  }
  return Yt;
}

grassmannQ::grassmannQ(int n1, int p1, int r1, 
          Rcpp::NumericMatrix Y1,int retraction1):
          manifold(n1,p1,r1,Y1,retraction1){}


double grassmannQ::evalHessian(const arma::mat &H, const arma::mat &Z){
  arma::mat projH=H;
  projH=projH-Y*Y.t()*projH;
  hessian_Z=projH+Z*Y.t()*xi;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}



void grassmannQ::vectorTrans(){
    descD=descD-Yt*Yt.t()*descD;
}

void grassmannQ::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

void grassmannQ::set_particle(){
  arma::mat y_temp=arma::randn(n,p);
  arma::mat Q,R;
  arma::qr_econ(Q,R,y_temp);
  Y=Q;
  arma::mat velocity_temp=arma::randn(n,p);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}

