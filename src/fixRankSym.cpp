#include "fixRankSym.h"


fixRankSym::fixRankSym(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){
               V=arma::mat(p,r,arma::fill::eye);
               Sigma=arma::vec(r,arma::fill::ones);
               Y=V*arma::diagmat(Sigma)*V.t();
//             arma::svd_econ(U,Sigma,V,Y);
//             U=U(arma::span::all,arma::span(0,r-1));
//             V=V(arma::span::all,arma::span(0,r-1));
//             Sigma=Sigma(arma::span(0,r-1));
           }


// V[V^T*gradF.Sym*V]V^T+  [(I-VV^T)*gradF.Sym*V]V^T+  V[V^T*gradF.Sym*(I-VV^T)]
// V[] M ]V^T+[ Up ]V^T+V[ Vp^T ]
void fixRankSym::evalGradient(arma::mat gradF,std::string method){
    arma::mat Rv;
    gradF=0.5*(gradF+gradF.t());
    Rv=gradF*V;
    M=V.t()*Rv;
    Vp=Rv-V*M;
    xi=V*M*V.t()+Vp*V.t()+V*Vp.t();
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
      M_desc=-M;
      Vp_desc=-Vp;
    }else if(method=="trustRegion"){
      xi_normal=gradF-xi;
    }
}


void fixRankSym::retrieve_steepest(){
  descD=-xi;
  M_desc=-M;
  Vp_desc=-Vp;
}

//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double fixRankSym::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat eucH_proj,Weingarten;
  //project euclidian Hessian onto tangent space
  arma::mat Rv,MM;
  arma::mat eucHSym=0.5*(eucH+eucH.t());
  //Ru=U.t()*eucH;
  Rv=eucHSym*V;
  MM=V.t()*Rv;
  Rv=(Rv-V*MM)*V.t();
  eucH_proj=V*MM*V.t()+Rv+Rv.t();
  //Weingarten Map
  arma::vec Sigma_inv=Sigma;
  int i;
  for(i=0;i<r;i++){
     Sigma_inv[i]=1/Sigma_inv[i];
  }
  arma::mat Y_plus=V*arma::diagmat(Sigma_inv)*V.t();
  Weingarten=xi_normal*Z.t()*Y_plus+Y_plus*Z.t()*xi_normal;
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}

//second argument unused
//rectract with step size
//first used to avoid computation repetition
arma::mat fixRankSym::retract(double stepSize, std::string method,bool first){
  arma::mat S,Vs,Us;
  arma::vec Sigma_s;
  if(first){
    arma::qr_econ(Qv,Rv,Vp_desc);
  }
  S=arma::mat(2*r,2*r,arma::fill::zeros);
  S(arma::span(0,r-1),arma::span(0,r-1))=arma::diagmat(Sigma)+M_desc*stepSize;
  S(arma::span(r,2*r-1),arma::span(0,r-1))=Rv*stepSize;
  S(arma::span(0,r-1),arma::span(r,2*r-1))=Rv.t()*stepSize; 
  //arma::eig_sym(Sigma_s,Vs,S);
  arma::svd_econ(Us,Sigma_s,Vs,S);
  Sigma_t=Sigma_s(arma::span(0,r-1));
  Vt=V*Vs(arma::span(0,r-1),arma::span(0,r-1))+Qv*Vs(arma::span(r,2*r-1),arma::span(0,r-1));
  Yt=Vt*arma::diagmat(Sigma_t)*Vt.t();
  return Yt;
}



void fixRankSym::acceptY(){
  Y=Yt;
  V=Vt;
  Sigma=Sigma_t;
}


//decompose xi_temp into VMV^T+ UpV^T+ VVp^T
 void fixRankSym::set_descD(arma::mat xi_temp){
  arma::mat Rv;//,Ru;
  Rv=xi_temp*V;
  M_desc=V.t()*Rv;
  Vp_desc=Rv-V*M_desc;
}



//vector transport of conjugate direction(update descent direction)
//from Y to Yt, evaluated before accept Y

void fixRankSym::vectorTrans(){
  arma::mat Av,Bv,Bu;
  Av=V.t()*Vt;  //r*r
  Bv=Vp_desc.t()*Vt;//r*r
  //Bu=Up_desc.t()*Vt;
  Vp_desc=V*M_desc*Av+Vp_desc*Av+V*Bv;
  //Vp_desc=V*M_desc.t()*Av+V*Bu+Vp_desc*Av;
  //M_desc=Av.t()*M_desc*Av+Bv.t()*Av+Av.t()*Bv;
  M_desc=Vt.t()*Vp_desc;
  Vp_desc=Vp_desc-Vt*M_desc;
  //Vp_desc=Vp_desc-Vt*Vt.t()*Vp_desc;
}


void fixRankSym::update_conjugateD(double eta){
  M_desc=M_desc*eta-M;
  //Up_desc=Up_desc*eta-Up;
  Vp_desc=Vp_desc*eta-Vp;
}

void fixRankSym::set_particle(){
//  arma::mat y_temp=arma::randn(n,p);
//  arma::mat Q,R;
//  arma::qr_econ(Q,R,y_temp);
//  Y=Q;
//  arma::mat velocity_temp=arma::randn(n,p);
//  //psedo gradient as velocity;
//  evalGradient(velocity_temp,"steepest");
//  conjugateD=xi;
}
