#include "fixRankPSD.h"


fixRankPSD::fixRankPSD(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){
               Y=arma::mat(n,r,arma::fill::eye);
           }



void fixRankPSD::evalGradient(arma::mat gradF,std::string method){
    //Omega=0 as gradF=dF*Y+dF^T*Y and Y^T*gradF=gradF^T*Y^T
    xi=gradF;
    if(method=="steepest"){
      eDescent=arma::dot(gradF,xi);
      descD=-xi;
    }
}


//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double fixRankPSD::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  arma::mat YtY=Y.t()*Y;
  arma::mat C=Y.t()*eucH;
  C=C.t()-C;
  hessian_Z=arma::syl(YtY,YtY,C);

  hessian_Z=eucH-Y*hessian_Z;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}

//second argument unused
//rectract with step size
//first used to avoid computation repetition
arma::mat fixRankPSD::retract(double stepSize, std::string method,
                  bool first){
  Yt=Y+stepSize*descD;
  return Yt;
}



//void fixRankPSD::acceptY(){
//  Y=Yt;
//}


// void fixRankPSD::set_descD(arma::mat xi_temp){
//descD=xi_temp;
//}
//



//vector transport of conjugate direction(update descent direction)
//from Y to Yt, evaluated before accept Y

void fixRankPSD::vectorTrans(){
    arma::mat A=Yt.t()*Yt;
    arma::mat C=Yt.t()*descD;
    C=C.t()-C;
    //solve the equation A*Omega+Omega*A+C=0
    arma::mat Omega=arma::syl(A,A,C);
    descD=descD-Yt*Omega;
}


void fixRankPSD::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

void fixRankPSD::set_particle(){

}
