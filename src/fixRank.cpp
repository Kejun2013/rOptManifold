#include "fixRank.h"


fixRank::fixRank(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){
               U=arma::mat(n,r,arma::fill::eye);
               V=arma::mat(p,r,arma::fill::eye);
               Sigma=arma::vec(r,arma::fill::ones);
               Y=U*arma::diagmat(Sigma)*V.t();
//             arma::svd_econ(U,Sigma,V,Y);
//             U=U(arma::span::all,arma::span(0,r-1));
//             V=V(arma::span::all,arma::span(0,r-1));
//             Sigma=Sigma(arma::span(0,r-1));
           }



void fixRank::evalGradient(arma::mat gradF,std::string method){
    arma::mat Ru,Rv;
    Ru=U.t()*gradF;
    Rv=gradF*V;
    M=U.t()*Rv;
    Up=Rv-U*M;
    Vp=(Ru-M*V.t()).t();
    xi=U*M*V.t()+Up*V.t()+U*Vp.t();
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
      M_desc=-M;
      Up_desc=-Up;
      Vp_desc=-Vp;
    }else if(method=="trustRegion"){
      //xi_normal=gradF-U*Ru-Rv*V.t()+U*M*V.t();
      xi_normal=gradF-xi;
    }else if(method=="particleSwarm"){
      descD=xi;
      M_desc=M;
      Up_desc=Up;
      Vp_desc=Vp;
    }
}


void fixRank::retrieve_steepest(){
  descD=-xi;
  M_desc=-M;
  Up_desc=-Up;
  Vp_desc=-Vp;
}

//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double fixRank::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat YU,eucH_proj,Weingarten;
  //project euclidian Hessian onto tangent space
  arma::mat Ru,Rv,MM;
  Ru=U.t()*eucH;
  Rv=eucH*V;
  MM=U.t()*Rv;
  eucH_proj=U*MM*V.t()+(Rv-U*MM)*V.t()+U*(Ru-MM*V.t());
  //Weingarten Map
  arma::vec Sigma_inv=Sigma;
  int i;
  for(i=0;i<r;i++){
     Sigma_inv[i]=1/Sigma_inv[i];
  }
  arma::mat Y_plus=V*arma::diagmat(Sigma_inv)*U.t();
  Weingarten=xi_normal*Z.t()*Y_plus.t()+Y_plus.t()*Z.t()*xi_normal;
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}

//second argument unused
//rectract with step size
//first used to avoid computation repetition
arma::mat fixRank::retract(double stepSize, std::string method,bool first){
  arma::mat S,Us,Vs;
  arma::vec Sigma_s;
  bool conv;
  Yt=Y;
  if(first){
    conv=arma::qr_econ(Qu,Ru,Up_desc);
    if(!conv) return Y;
    conv=arma::qr_econ(Qv,Rv,Vp_desc);
    if(!conv) return Y;
  }
  S=arma::mat(2*r,2*r,arma::fill::zeros);
  S(arma::span(0,r-1),arma::span(0,r-1))=arma::diagmat(Sigma)+M_desc*stepSize;
  S(arma::span(r,2*r-1),arma::span(0,r-1))=Ru*stepSize;
  S(arma::span(0,r-1),arma::span(r,2*r-1))=Rv.t()*stepSize; 
  conv=arma::svd_econ(Us,Sigma_s,Vs,S);
  if(!conv) return Y;
  Sigma_t=Sigma_s(arma::span(0,r-1));
  Ut=U*Us(arma::span(0,r-1),arma::span(0,r-1))+Qu*Us(arma::span(r,2*r-1),arma::span(0,r-1));
  Vt=V*Vs(arma::span(0,r-1),arma::span(0,r-1))+Qv*Vs(arma::span(r,2*r-1),arma::span(0,r-1));
  
  Yt=Ut*arma::diagmat(Sigma_t)*Vt.t();
  return Yt;
}



void fixRank::acceptY(){
  Y=Yt;
  U=Ut;
  V=Vt;
  Sigma=Sigma_t;
}


 void fixRank::set_descD(arma::mat xi_temp){
  arma::mat Rv,Ru;
  Rv=xi_temp*V;
  Ru=U.t()*xi_temp;
  M_desc=Ru*V;
  Up_desc=Rv-U*M_desc;
  Vp_desc=(Ru-M_desc*V.t()).t();
}



//vector transport of conjugate direction(update descent direction)
//from Y to Yt, evaluated before accept Y

void fixRank::vectorTrans(){
  arma::mat Av,Au,Bv,Bu;
  Av=V.t()*Vt;
  Au=U.t()*Ut;
  Bv=Vp_desc.t()*Vt;
  Bu=Up_desc.t()*Ut;
  Up_desc=U*M_desc*Av+Up_desc*Av+U*Bv;
  Vp_desc=V*M_desc.t()*Au+V*Bu+Vp_desc*Au;
  M_desc=Au.t()*M_desc*Av+Bu.t()*Av+Au.t()*Bv;
  Up_desc=Up_desc-Ut*Ut.t()*Up_desc;
  Vp_desc=Vp_desc-Vt*Vt.t()*Vp_desc;
}


void fixRank::update_conjugateD(double eta){
  M_desc=M_desc*eta-M;
  Up_desc=Up_desc*eta-Up;
  Vp_desc=Vp_desc*eta-Vp;
}

void fixRank::set_particle(){ 
  U=arma::randn(n,r);
  V=arma::randn(p,r);
  arma::mat Q,R;
  arma::qr_econ(Q,R,U);
  U=Q;
  arma::qr_econ(Q,R,V);
  V=Q;
  Sigma=arma::vec(r,arma::fill::ones);
  Y=U*arma::diagmat(Sigma)*V.t();
  //////////////////
  arma::mat velocity_temp=arma::randn(n,p);
  evalGradient(velocity_temp,"steepest");
}
