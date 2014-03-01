#include <omp.h>
#include "rcpparma_hello_world.h"

using namespace Rcpp ;



SEXP rcpparma_hello_world(SEXP n1, SEXP p1, SEXP f1, SEXP f2, SEXP f3, 
                          SEXP iterMax1,SEXP retract1){
  int retract=as< int>(retract1);//0 expoential 1 QR 2 Cayley
	int n,p,iter=0,iterMax=as< int>(iterMax1),i;
  Function obej(f1);
  Function grad(f2);
  Function expm(f3);
  n=as< int>(n1);
  p=as< int>(p1);
  //Function f(f1);
	arma::mat Y = arma::eye<arma::mat>(n, p);
  arma::mat Delta(n,p), Yt(n,p),Yt_best(n,p),gradF;
  arma::mat QR,Q,R,K(2*p,2*p,arma::fill::zeros);
	NumericMatrix gradF2;
  double step=0,objValue=as< double>(obej(Y)), objValue_temp;
  //begin iteration
  while(iter<iterMax){
    //gradient of objective funtion
    iter++;
    gradF2=grad(Y);
    gradF=arma::mat(gradF2.begin(),n,p,false);
    //gradient on the stiefel manifold
    Delta=gradF-Y*gradF.t()*Y;
    //  //prepare for matrix expoential
    if(retract==0){
      QR=arma::eye<arma::mat>(n,n)-Y*Y.t();
      QR=QR*Delta;
      arma::qr_econ(Q,R,QR);
      K(arma::span(0,p-1),arma::span(0,p-1))=Y.t()*Delta;
      K(arma::span(0,p-1),arma::span(p,2*p-1))=-R.t();
      K(arma::span(p,2*p-1),arma::span(0,p-1))=R;
    }
    //prepare for Cayley 
    if(retract==2){
      Q=arma::mat(n,2*p,arma::fill::zeros);
      R=Q;
      Q(arma::span::all,arma::span(0,p-1))=gradF;
      Q(arma::span::all,arma::span(p,2*p-1))=Y;
      R(arma::span::all,arma::span(0,p-1))=Y;
      R(arma::span::all,arma::span(p,2*p-1))=-gradF;
    }
    //step size selection
    step=10;Yt_best=Y;
    for(i=0;i<20;i++){
      step=step*0.66;
      if(retract==1){
        Yt=Y-step*Delta;
        arma::qr_econ(Q,R,Yt);
        Yt=Q;
      }else if(retract==0){
        gradF2=expm(K*step);
        gradF=arma::mat(gradF2.begin(),2*p,2*p,false);
        Yt=Y*gradF(arma::span(0,p-1),arma::span(0,p-1))+Q*gradF(arma::span(p,2*p-1),arma::span(0,p-1));
      }else{
        Yt=arma::eye<arma::mat>(2*p,2*p)+step/2*R.t()*Q;
        Yt=Y-step*Q*Yt.i()*R.t()*Y;
      }
      objValue_temp=as< double>(obej(Yt));
      if(objValue_temp<objValue){
        objValue=objValue_temp;
        Yt_best=Yt;
      }
    }
    Y=Yt_best;
  }
  return List::create(Named("optY")=Y,Named("optValue")=objValue);
}

