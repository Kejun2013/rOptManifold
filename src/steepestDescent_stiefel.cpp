//#include <omp.h>
#include "steepestDescent_stiefel.h"

using namespace Rcpp;
using namespace std;

double qr_stiefel(List &YList, vector<arma::mat> &Y, const int k, const int prodK,
              arma::mat Delta,
              Function obej, double &objValue){
    double  step=5,stepsize_best=0, objValue_best=objValue,objValue_temp;
    arma::mat Q,R,Yt,Yt_best=Y[k];
    int i;
    for(i=0;i<40;i++){//c
      step=step*0.66;
      Yt=Y[k]-step*Delta;
      arma::qr_econ(Q,R,Yt);
      if(R(0,0)<0){
        Q=-Q;
      }
      YList[k]=Q;
      if(prodK>1){
        objValue_temp=as< double>(obej(YList));
      }else{
        objValue_temp=as< double>(obej(YList[0]));
      }
      if(objValue_temp<objValue_best){
        objValue_best=objValue_temp;
        Yt_best=Q;
        stepsize_best=step;
      }
    }
    objValue=objValue_best;
    Y[k]=Yt_best;
    YList[k]=Yt_best;
    return stepsize_best;
}

double Cayley_stiefel(List &YList, vector<arma::mat> &Y, const int k,
                   const int prodK,const int n, const int p, 
                  arma::mat gradF, Function obej, double &objValue){
      arma::mat Q,R,Yt,Yt_best=Y[k];
      Q=arma::mat(n,2*p,arma::fill::zeros);
      R=Q;
      Q(arma::span::all,arma::span(0,p-1))=gradF;
      Q(arma::span::all,arma::span(p,2*p-1))=Y[k];
      R(arma::span::all,arma::span(0,p-1))=Y[k];
      R(arma::span::all,arma::span(p,2*p-1))=-gradF;
      double  step=5,stepsize_best=0, objValue_best=objValue,objValue_temp;
    int i;
    for(i=0;i<40;i++){//c
      step=step*0.66;
      Yt=arma::eye<arma::mat>(2*p,2*p)+step/2*R.t()*Q;
      Yt=Y[k]-step*Q*Yt.i()*R.t()*Y[k];
      YList[k]=Yt;
      if(prodK>1){
        objValue_temp=as< double>(obej(YList));
      }else{
        objValue_temp=as< double>(obej(YList[0]));
      }
      if(objValue_temp<objValue_best){
        objValue_best=objValue_temp;
        Yt_best=Yt;
        stepsize_best=step;
      }
    }
    objValue=objValue_best;
    Y[k]=Yt_best;
    YList[k]=Yt_best;
    return stepsize_best;
}

double expm_stiefel(List &YList, vector<arma::mat> &Y, 
              const int k, const int prodK, const int n, const int p,
              arma::mat Delta, Function obej, Function expm, double &objValue){
    arma::mat Q,R,QR, gradF,Yt,Yt_best=Y[k];
    NumericMatrix gradF2;
    QR=arma::eye<arma::mat>(n,n)-Y[k]*Y[k].t();
    QR=QR*Delta;
    arma::qr_econ(Q,R,QR);
    arma::mat K(2*p,2*p,arma::fill::zeros);
    K(arma::span(0,p-1),arma::span(0,p-1))=Y[k].t()*Delta;
    K(arma::span(0,p-1),arma::span(p,2*p-1))=-R.t();
    K(arma::span(p,2*p-1),arma::span(0,p-1))=R;
    double  step=5,stepsize_best=0, objValue_best=objValue,objValue_temp;
    int i;
    for(i=0;i<40;i++){//c
      step=step*0.66;
      gradF2=expm(K*step);
      gradF=arma::mat(gradF2.begin(),2*p,2*p,false);
      Yt=Y[k]*gradF(arma::span(0,p-1),arma::span(0,p-1))+Q*gradF(arma::span(p,2*p-1),arma::span(0,p-1));
      YList[k]=Yt;
      if(prodK>1){
        objValue_temp=as< double>(obej(YList));
      }else{
        objValue_temp=as< double>(obej(YList[0]));
      }
      if(objValue_temp<objValue_best){
        objValue_best=objValue_temp;
        Yt_best=Yt;
        stepsize_best=step;
      }
    }
    objValue=objValue_best;
    Y[k]=Yt_best;
    YList[k]=Yt_best;
    return stepsize_best;
}

SEXP steepestDescent_stiefel(SEXP n1, SEXP p1, SEXP f1, SEXP f2, SEXP f3, 
                          SEXP iterMax1,SEXP tol1,SEXP retract1){
  int retract=as< int>(retract1);//0 expoential 1 QR 2 Cayley
	int iter=0,iterMax=as< int>(iterMax1),k;
  double tol=as< double>(tol1),objValue;
  Function obej(f1);
  Function grad(f2);
  Function expm(f3);
  NumericVector n(n1), p(p1);
  int prodK=n.size();// size of product manifold
	vector< arma::mat> Y(prodK);
  List YList(prodK);
  for(k=0;k<prodK;k++){
     Y[k]= arma::eye<arma::mat>(n[k], p[k]);
     YList[k]=Y[k];
  }
  arma::mat Delta,gradF; //Delta: steepest tangent
  //arma::mat QR,Q,R,K(2*p,2*p,arma::fill::zeros);
	NumericMatrix gradF2;
  if(prodK>1){
    objValue=as< double>(obej(YList));
  }else{
    objValue=as< double>(obej(YList[0]));
  }
  //begin iteration
  bool flag=true;
  vector<double> stepsize(prodK);
  while(iter<iterMax && flag){
    //gradient of objective funtion
    iter++;
    for(k=0;k<prodK;k++){
        if(prodK>1){
          gradF2=grad(YList,k+1);
        }else{
          gradF2=grad(YList[0]);
        }
        gradF=arma::mat(gradF2.begin(),n[k],p[k],false);
        //gradient on the stiefel manifold
        if(retract==1){
          Delta=gradF-Y[k]*gradF.t()*Y[k];
          stepsize[k]=qr_stiefel(YList,Y,k,prodK,Delta,obej,objValue);
        }else if(retract==2){
          stepsize[k]=Cayley_stiefel(YList, Y, k,prodK, n[k], p[k], 
                           gradF,  obej, objValue);
        }else{
          Delta=gradF-Y[k]*gradF.t()*Y[k];
          stepsize[k]=expm_stiefel(YList, Y,k, prodK, n[k],p[k],
                  Delta, obej, expm, objValue);
        }
    }
    if(tol>(*max_element(stepsize.begin(),stepsize.end()))) flag=false;
  }
  return List::create(Named("optY")=YList,Named("optValue")=objValue);
}//end of function



