//#include <omp.h>
#include "stiefel.h"
#include "grassmanQ.h"

// YList is a (list of) of matrix(ces) of dimensional n1*p1, initial values
//f1 is the objective function; f2 is the gradient function; f3 is the hessian function
// iterMax: maximum number of iterations
// tol: tolerance of accuracy
//retraction: retraction method

RcppExport SEXP  trustRegion(SEXP YList1, SEXP n1, SEXP p1, SEXP r1,
                      SEXP mtype1,SEXP retraction1,
                      SEXP f1, SEXP f2,SEXP f3,
                      SEXP control1){
BEGIN_RCPP
  //Initialization of functions and control parameters
  Function obej(f1);
  Function grad(f2);
  Function hessian(f3);
  List control(control1);
  IntegerVector retraction(retraction1);
  int iterMax=as< int>(control["iterMax"]);
  double tol=as< double>(control["tol"]);
  double Delta0=as< double>(control["delta"]);
  double DeltaMax=as< double>(control["deltaMax"]);
  double theta=as< double>(control["theta"]);
  double kappa=as< double>(control["kappa"]);
  double rhoMin=as< double>(control["rhoMin"]);
  
  // Initialization of Data points
  IntegerVector n(n1),p(p1),r(r1);
  CharacterVector mtype(mtype1);
  int prodK=n.size();// size of product manifold
  vector< manifold*>  manifoldY;
  vector< double> Delta(prodK,Delta0);
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
     else if(typeTemp=="grassmanQ"){
       manifoldY.push_back(new grassmanQ(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="fixedRank"){
       
     }else if(typeTemp=="grassmanS")
  }
    
  //define other varibles
  int iter=0,iterSub; // outer loop, inner loop control
  bool flag=true,bdryTouched;
  double objValue,objValue_temp,objDesc,m_temp,mDesc;
  double alpha,beta,qA,qB,qC,tau,rrNorm0,rrNorm,
          etaNorm,etaNorm_new,hessQterm;
  //value of objective function
  //eDescent: expected descent amount
  //largest obj descent amount
//  Function expm(f3);
  arma::mat gradF,r,d,eta;  //gradient in ambient space
  
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
         manifoldY[k]->evalGradient(gradF,"trustRegion");
         
        r=manifoldY[k]->getGradient();
        d=-r; //little delta
        eta=arma::mat(n[k],p[k],arma::fill::zeros);
        rrNorm=manifoldY[k]->metric(r,r);  //<r,r>_Y
        //rrNorm0: controlling stopping condition for the subproblem
        rrNorm0=rrNorm;
        rrNorm0=rrNorm0*min(pow(rrNorm0,theta),kappa);
        etaNorm=0;
        //begin of inner loop: sub-problem
        iterSub=0;
        bdryTouched=false;
        while(iterSub<1000){
          iterSub++;
          if(prodK>1){
             hessianF=as< arma::mat>(hessian(YList,k+1,d));
          }else{
            hessianF=as< arma::mat>(hessian(YList[0],d));
          }
          //hessQterm=<delta_j,H_k*delta_j>
          hessQterm=manifoldY[k]->evalHessian(hessianF,d);
          if(hessQterm<=0.0001){
            //solve a quadratic function to find the point touching boundary
            qA=manifoldY[k]->metric(d,d);
            qB=manifoldY[k]->metric(eta,d);
            qC=-Delta[k]*Delta[k]+etaNorm;
            tau=sqrt(qB*qB-4*qA*qC);
            tau=(tau-qB)/(2*qA);
            bdryTouched=true;
            break;
          }
          alpha=rrNorm/hessQterm;
          eta=eta+alpha*d;
          etaNorm_new=manifoldY[k]->metric(eta,eta);
          if(etaNorm_new>0.9999*Delta){
            qA=manifoldY[k]->metric(d,d);
            qB=manifoldY[k]->metric(eta,d);
            qC=-Delta[k]*Delta[k]+etaNorm;
            tau=sqrt(qB*qB-4*qA*qC);
            tau=(tau-qB)/(2*qA);
            bdryTouched=true;
            break;
          }
          etaNorm=etaNorm_new;
          r=r+alpha*(manifoldY[k]->get_hessianZ());
          beta=(manifoldY[k]->metric(r,r))/rrNorm;
          rrNorm=beta*rrNorm;
          if(rrNorm<rrNorm0)){
            break;
          }
          d=-r+beta*d;
        }
        //end of inner loop:: subproblem
        m_temp=manifoldY[k]->secondOrderApprox(objValue,eta);
        manifoldY[k]->set_Gradient(eta);
        YList[k]=manifoldY[k]->retract(1);
        if(prodK>1){
          objValue_temp=as< double>(obej(YList));
        }else{
          objValue_temp=as< double>(obej(YList[0]));
        }
        objDesc=objValue-objValue_temp;
        mDesc=objValue-m_temp;
        rho=objDesc/mDesc;
        if(rho<0.25){
          Delta=0.25*Delta;
        }else if(rho>0.75 && bdryTouched){
          Delta=min(2*Delta,DeltaMax);
        }
        if(rho>rhoMin){
          manifoldY[k]->acceptY();
          objValue=objValue_temp;
        }else{
          YList[k]=manifoldY[k]->get_Y();
        }
    }// iteration over product component
    if(tol> (*max_element(Delta.begin(),Delta.end()))) flag=false;
  }// outer iteration
  return List::create(Named("optY")=YList,
              Named("optValue")=objValue,
              Named("NumIter")=iter);
END_RCPP
}//end of function
