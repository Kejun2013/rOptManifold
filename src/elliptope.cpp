#include "elliptope.h"


elliptope::elliptope(int n1, int p1, int r1, 
           Rcpp::NumericMatrix Y1,int retraction1):
           manifold(n1,p1,r1,Y1,retraction1){
//               Y=arma::mat(n,r,arma::fill::eye);
//               int i;
//               for(i=r;i<n;i++){
//                 Y(i,1)=1;
//               }
           }



void elliptope::evalGradient(arma::mat gradF,std::string method){
    //Omega=0 as gradF=dF*Y+dF^T*Y and Y^T*gradF=gradF^T*Y^T
    arma::mat temp=gradF%Y;
    //alpha_i=tr(gradF^T*A_i*Y)
    alphas=arma::sum(temp,1);
    xi=gradF-(alphas*arma::mat(1,r,arma::fill::ones))%Y;
    if(method=="steepest"){
      descD=-xi;
      eDescent=arma::dot(xi,gradF);
    }else{
      eGrad=gradF;//storage for hessian
    }
}



//evaluate Hessian operator on Z
//eucH is the ordinary Hessian matrix on Z in the ambient space
//return <Z, Hessian*Z>_Y
double elliptope::evalHessian(const arma::mat & eucH,const arma::mat & Z){
  
  hessian_Z=eucH-(alphas*arma::mat(1,r,arma::fill::ones))%Z;
  arma::colvec DaZ=arma::sum(Z%Y,1);
  DaZ=(-2.0)*DaZ%DaZ;
  DaZ=DaZ+arma::sum(eucH%Y+eGrad%Z,1);
  hessian_Z=hessian_Z-(alphas*arma::mat(1,r,arma::fill::ones))%Y;
  //projection onto horizontal space
  arma::mat YtY=Y.t()*Y;
  arma::mat C=Y.t()*hessian_Z;
  C=C.t()-C;
  arma::mat Omega=arma::syl(YtY,YtY,C);
  arma::colvec betas=arma::sum(hessian_Z%Y,1);
  hessian_Z=hessian_Z-Y*Omega-(betas*arma::mat(1,r,arma::fill::ones))%Y;
  Z_hessian_Z=arma::dot(hessian_Z,Z);
  return Z_hessian_Z;
}

//second argument unused
//rectract with step size
//first used to avoid computation repetition
arma::mat elliptope::retract(double stepSize, std::string method,
                  bool first){
  if(retraction==0){//exponential
    if(first){
      retract_a=arma::sum(descD%descD,1);
      retract_a=arma::sqrt(retract_a);
    }
    arma::mat mRepeat=arma::mat(1,r,arma::fill::ones);
    Yt=(arma::cos(stepSize*retract_a)*mRepeat)%Y+
          (1/stepSize*arma::pow(retract_a,-1))%
          (arma::sin(stepSize*retract_a)*mRepeat)%xi;
  }else{//projection
    Yt=Y+stepSize*descD;
    retract_a=arma::sum(Yt%Yt,1);
    retract_a=arma::pow(retract_a,-0.5);
    Yt=(retract_a*arma::mat(1,r,arma::fill::ones))%Yt;
  }
  return Yt;
}





//vector transport of conjugate direction(update descent direction)
//from Y to Yt, evaluated before accept Y

void elliptope::vectorTrans(){
  arma::mat YtY=Yt.t()*Yt;
  arma::mat C=Yt.t()*descD;
  C=C.t()-C;
  arma::mat Omega=arma::syl(YtY,YtY,C);
  arma::colvec betas=arma::sum(descD%Yt,1);
  descD=descD-Yt*Omega-(betas*arma::mat(1,r,arma::fill::ones))%Yt;
}


void elliptope::update_conjugateD(double eta){
  descD=eta*descD-xi;
}

void elliptope::set_particle(){
      Y=arma::randn(n,r);
  arma::mat velocity_temp=arma::randn(n,r);
  //psuedo gradient as velocity;
  evalGradient(velocity_temp,"particleSwarm");
}
