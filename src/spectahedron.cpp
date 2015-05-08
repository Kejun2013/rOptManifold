#include "spectahedron.h"


spectahedron::spectahedron(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){
               Y=arma::mat(n,r,arma::fill::eye);
               Y=Y/sqrt(r);
           }



void spectahedron::evalGradient(arma::mat gradF,std::string method){
    //Omega=0 as gradF=dF*Y+dF^T*Y and Y^T*gradF=gradF^T*Y^T
    alpha=arma::dot(gradF.t(),Y);
    xi=gradF-alpha*Y;
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
    }else if(method=="trustRegion"){
      eGrad=gradF;
    }else if(method=="particleSwarm"){
      descD=xi;
    }
}



//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double spectahedron::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  hessian_Z=eucH-alpha*Z-
              (arma::dot(eucH,Y)+
              arma::dot(eGrad,Z)-pow(arma::dot(Z,Y),2))*Y;
  arma::mat YtY=Y.t()*Y;
  arma::mat C=Y.t()*hessian_Z;
  C=C.t()-C;
  arma::mat Omega=arma::syl(YtY,YtY,C);
  double beta=arma::dot(Y,hessian_Z);
  hessian_Z=hessian_Z-Y*Omega-beta*Y;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}

//second argument unused
//rectract with step size
//first used to avoid computation repetition
arma::mat spectahedron::retract(double stepSize, std::string method,
                  bool first){
  if(retraction==0){//exponential
    if(first){
      retract_a=arma::norm(descD,"fro");
    }
    Yt=cos(stepSize*retract_a)*Y+sin(stepSize*retract_a)*descD/retract_a;
    //Yt=Yt/arma::norm(Yt,"fro");
  }else{//projection
    Yt=Y+stepSize*descD;
    retract_a=arma::norm(Yt,"fro");
    Yt=Yt/retract_a;
  }
  return Yt;
}





//vector transport of conjugate direction(update descent direction)
//from Y to Yt, evaluated before accept Y

void spectahedron::vectorTrans(){
    arma::mat A=Yt.t()*Yt;
    arma::mat C=Yt.t()*descD;
    C=C.t()-C;
    //solve the equation A*Omega+Omega*A+C=0
    arma::mat Omega=arma::syl(A,A,C);
    double beta=arma::dot(descD,Yt);
    descD=descD-Yt*Omega-beta*Y;
}


void spectahedron::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

void spectahedron::set_particle(){
  Y=arma::randn(n,r);
  double normTemp=arma::norm(Y,"fro");
  Y=Y/normTemp;
  //psuedo gradient as velocity;
   arma::mat velocity_temp=arma::randn(n,r);
  evalGradient(velocity_temp,"particleSwarm");
}
