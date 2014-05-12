//#include <omp.h>
#include "stiefel.h"
#include "grassmannQ.h"

// YList is a (list of) of matrix(ces) of dimensional n1*p1, initial values
//f1 is objective function; f2 is gradient function
// iterMax: maximum number of iterations
// tol: tolerance of accuracy
//retract: retraction method

RcppExport SEXP  steepestDescent(SEXP YList1, SEXP n1, SEXP p1, SEXP r1,
                      SEXP mtype1,SEXP retraction1,
                      SEXP f1, SEXP f2,
                      SEXP control1){
BEGIN_RCPP
  //Initialization of functions and control parameters
  Function obej(f1);
  Function grad(f2);
  List control(control1);
  IntegerVector retraction(retraction1);
  int iterMax=as< int>(control["iterMax"]);
  double tol=as< double>(control["tol"]);
  double sigma=as< double>(control["sigma"]);
  double beta=as< double>(control["beta"]);
  double alpha=as< double>(control["alpha"]);

  
  // Initialization of Data points
  IntegerVector n(n1),p(p1),r(r1);
  CharacterVector mtype(mtype1);
  int prodK=n.size();// size of product manifold
  vector< manifold*>  manifoldY;
  int k;
  List YList(YList1);
  for(k=0;k<prodK;k++){
    SEXP yTemp2=YList[k];
     NumericMatrix yTemp(yTemp2);
     std::string typeTemp=as< std::string>(mtype[k]);
     if(typeTemp=="stiefel"){
        manifoldY.push_back(new stiefel(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }
     else if(typeTemp=="grassmannQ"){
       manifoldY.push_back(new grassmannQ(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="fixedRank"){
       
     }
  }
    
  //define other varibles
	int iter=0; // outer loop, inner loop control
  bool flag=true;
  double objValue,objValue_temp,eDescent,objDesc,stepsize;
  //value of objective function
  //eDescent: expected descent amount
  //largest obj descent amount
//  Function expm(f3);
  arma::mat gradF;  //gradient in ambient space
  
  if(prodK>1){
    objValue=as< double>(obej(YList));
  }else{
    objValue=as< double>(obej(YList[0]));
  }
  
  //begin iteration
  while(iter<iterMax && flag){
    //gradient of objective funtion
    iter++;
    objDesc=-1;
    for(k=0;k<prodK;k++){
        if(prodK>1){
          gradF=as< arma::mat>(grad(YList,k+1));
        }else{
          gradF=as< arma::mat>(grad(YList[0]));
        }
        //gradient on the stiefel manifold
        manifoldY[k]->evalGradient(gradF,"steepest");
        stepsize=alpha/beta;
        eDescent=sigma/beta*(manifoldY[k]->get_eDescent());
        do{//c
          stepsize=stepsize*beta;
          eDescent=eDescent*beta;
          YList[k]=manifoldY[k]->retract(stepsize,"steepest");
          if(prodK>1){
            objValue_temp=as< double>(obej(YList));
          }else{
            objValue_temp=as< double>(obej(YList[0]));
          }
        }while((objValue-objValue_temp)<eDescent);
        //step size iteration
        //if a stepsize is accepted, update current location
        manifoldY[k]->acceptY();
        objDesc=objValue-objValue_temp;
        objValue=objValue_temp;
    }// iteration over product component
    if(tol>objDesc) flag=false;
  }// outer iteration
  return List::create(Named("optY")=YList,
              Named("optValue")=objValue,
              Named("NumIter")=iter);
END_RCPP
}//end of function




///////////////////////////////////////////////////////////////////////


//double expm_stiefel(List &YList, vector<arma::mat> &Y, 
//              const int k, const int prodK, const int n, const int p,
//              arma::mat Delta, Function obej, Function expm, double &objValue){
//    arma::mat Q,R,QR, gradF,Yt,Yt_best=Y[k];
//    NumericMatrix gradF2;
//    QR=arma::eye<arma::mat>(n,n)-Y[k]*Y[k].t();
//    QR=QR*Delta;
//    arma::qr_econ(Q,R,QR);
//    arma::mat K(2*p,2*p,arma::fill::zeros);
//    K(arma::span(0,p-1),arma::span(0,p-1))=Y[k].t()*Delta;
//    K(arma::span(0,p-1),arma::span(p,2*p-1))=-R.t();
//    K(arma::span(p,2*p-1),arma::span(0,p-1))=R;
//    double  step=5,stepsize_best=0, objValue_best=objValue,objValue_temp;
//    int i;
//    for(i=0;i<40;i++){//c
//      step=step*0.66;
//      gradF2=expm(K*step);
//      gradF=arma::mat(gradF2.begin(),2*p,2*p,false);
//      Yt=Y[k]*gradF(arma::span(0,p-1),arma::span(0,p-1))+Q*gradF(arma::span(p,2*p-1),arma::span(0,p-1));
//      YList[k]=Yt;
//      if(prodK>1){
//        objValue_temp=as< double>(obej(YList));
//      }else{
//        objValue_temp=as< double>(obej(YList[0]));
//      }
//      if(objValue_temp<objValue_best){
//        objValue_best=objValue_temp;
//        Yt_best=Yt;
//        stepsize_best=step;
//      }
//    }
//    objValue=objValue_best;
//    Y[k]=Yt_best;
//    YList[k]=Yt_best;
//    return stepsize_best;
//}


