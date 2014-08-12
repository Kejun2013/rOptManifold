#include "grassmannSub.h"

void grassmannSub::evalGradient(arma::mat gradF, std::string method){
  xi=Y*gradF-gradF*Y;
  xi=Y*xi-xi*Y;
  eDescent=arma::dot(gradF,xi);
  if(method=="steepest"){
    descD=-xi;
  }else if(method=="particleSwarm"){
    descD=xi;
  }else if(method=="trustRegion"){
    xi_normal=gradF-xi;
  }
}



arma::mat grassmannSub::retract(double stepSize, std::string method,bool first, Function expm){
  if(first){
      xiP=descD*Y-Y*descD;
  }
  if(retraction==1){ //QR
    arma::mat temp,temp_Q,temp_R;
    temp.eye(n,n);
    temp=temp+stepSize*xiP;
    arma::qr_econ(temp_Q,temp_R,temp);
    if(temp_R(0,0)<0){
      temp_Q=-temp_Q;
    }
    Yt=temp_Q*Y*temp_Q.t();
  }else if(retraction==2){//Cayley
    arma::mat temp1,temp2,temp3;
    temp1.eye(n,n);
    temp2.eye(n,n);
    temp1=2*temp1-xiP*stepSize;
    temp2=2*temp2+xiP*stepSize;
    temp3=temp2*temp1.i();
    temp2=temp1*temp2.i();
    Yt=temp3*Y*temp2;
  }
  return Yt;
}

grassmannSub::grassmannSub(int n1, int p1, int r1, 
          Rcpp::NumericMatrix Y1,int retraction1):
          manifold(n1,p1,r1,Y1,retraction1){}


double grassmannSub::evalHessian(const arma::mat &H, const arma::mat &Z){
  arma::mat eucH_proj,Weingarten;
  eucH_proj=Y*H-H*Y;
  eucH_proj=Y*eucH_proj-eucH_proj*Y;
  //
  Weingarten=xi_normal*Z-Z*xi_normal;
  Weingarten=Y*Weingarten-Weingarten*Y;
  Weingarten=-Weingarten;
  //
  hessian_Z=eucH_proj+Weingarten;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}


// ad^2_Yt(descD)
void grassmannSub::vectorTrans(){
    descD=Yt*descD-descD*Yt;
    descD=Yt*descD-descD*Yt;
}


void grassmannSub::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

void grassmannSub::set_particle(){
  arma::mat y_temp=arma::randn(n,r);
  arma::mat Q,R;
  arma::qr_econ(Q,R,y_temp);
  Y=Q*Q.t();
  arma::mat velocity_temp=arma::randn(n,n);
  //psedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
