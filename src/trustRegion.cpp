#include "stiefel.h"
#include "grassmannQ.h"
#include "grassmannSub.h"
#include "fixRank.h"
#include "fixRankPSD.h"
#include "fixRankSym.h"
#include "spectahedron.h"
#include "elliptope.h"
#include "sphere.h"
#include "oblique.h"

// YList is a (list of) of matrix(ces) of dimensional n1*p1, initial values
//f1 is the objective function; f2 is the gradient function; f3 is the hessian function
//mtype: type of manifold
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
  int iterSubMax=as< int>(control["iterSubMax"]);
  double tol=as< double>(control["tol"]);
  double Delta0=as< double>(control["Delta0"]);
  double DeltaMax=as< double>(control["DeltaMax"]);
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
     else if(typeTemp=="grassmannQ"){
       manifoldY.push_back(new grassmannQ(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="grassmannSub"){
       manifoldY.push_back(new grassmannSub(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="fixedRank"){
       manifoldY.push_back(new fixRank(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="fixedRankPSD"){
       manifoldY.push_back(new fixRankPSD(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="elliptope"){
       manifoldY.push_back(new elliptope(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="spectahedron"){
       manifoldY.push_back(new spectahedron(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="sphere"){
       manifoldY.push_back(new sphere(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="fixedRankSym"){
       manifoldY.push_back(new fixRankSym(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="oblique"){
       manifoldY.push_back(new oblique(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }       
  }
    
  //define other varibles
  int iter=0,iterSub; // outer loop, inner loop control
  bool flag=true,bdryTouched;
  
  //value of objective function, m function, related descent amount
  double objValue,objValue_temp,objValue_outer,objDesc,m_temp,mDesc;//,DeltaL;
  double alpha,beta,rho, qA,qB,qC,tau,rrNorm0,rrNorm,rrNormMax,
          etaNorm,etaNorm_new,hessQterm;
  //gradient, hessian of objective F in ambient space
  arma::mat gradF,hessianF, rr,dd,eta,eta_old,HZz;  ////////////////////////
  
  if(prodK>1){
    objValue=as< double>(obej(YList));
  }else{
    objValue=as< double>(obej(YList[0]));
  }
  
  //begin iteration
  while(iter<iterMax && flag){
    objValue_outer=objValue;
    //gradient of objective funtion
    iter++;
    rrNormMax=0;
    for(k=0;k<prodK;k++){
        if(prodK>1){
          gradF=as< arma::mat>(grad(YList,k+1));
        }else{
          gradF=as< arma::mat>(grad(YList[0]));
        }
         manifoldY[k]->evalGradient(gradF,"trustRegion");
         
        rr=manifoldY[k]->get_Gradient();
        dd=-rr; //little delta, n-by-p
        eta=arma::mat(dd.n_rows,dd.n_cols,arma::fill::zeros);  //optimizer of sub-problem
        rrNorm=manifoldY[k]->gradMetric();  //<r,r>_Y
        //rrNorm0: controlling stopping condition for the subproblem
        rrNorm0=rrNorm;
        rrNormMax=max(rrNorm,rrNormMax);
        rrNorm0=rrNorm0*min(pow(rrNorm0,theta),kappa);
        etaNorm=0;
        //begin of inner loop: sub-problem
        iterSub=0;
        bdryTouched=false;  //boundary of trust region touched?
        while(iterSub<iterSubMax){//////////////////
          iterSub++;
          if(prodK>1){
             hessianF=as< arma::mat>(hessian(YList,dd,k+1));
          }else{
            hessianF=as< arma::mat>(hessian(YList[0],dd));
          }
          hessQterm=manifoldY[k]->evalHessian(hessianF,dd);//<Z,H*Z>_Y, the quadratic term
          if(hessQterm<=0.0001){
            //delta is the direction of negative curvature
            //solve a quadratic equation to find the point touching boundary
            //qA*tau^2+qB*tau+qC=0 with unkown tau
            qA=manifoldY[k]->metric(dd,dd);
            qB=2*(manifoldY[k]->metric(eta,dd));
            qC=-Delta[k]*Delta[k]+etaNorm;
            tau=sqrt(qB*qB-4*qA*qC);
            tau=(tau-qB)/(2*qA);
            bdryTouched=true;
            eta=eta+tau*dd;
            break;
          }
          alpha=rrNorm/hessQterm;
          eta_old=eta;
          eta=eta+alpha*dd;
          etaNorm_new=manifoldY[k]->metric(eta,eta);
          if(etaNorm_new>0.9999*Delta[k]*Delta[k]){
            //if the optimizer is out of trust region
            qA=manifoldY[k]->metric(dd,dd);
            qB=2*(manifoldY[k]->metric(eta_old,dd));
            qC=-Delta[k]*Delta[k]+etaNorm;
            tau=sqrt(qB*qB-4*qA*qC);
            tau=(tau-qB)/(2*qA);
            eta=eta_old+tau*dd;
            bdryTouched=true;
            break;
          }
           etaNorm=etaNorm_new;
          rr=rr+alpha*(manifoldY[k]->get_hessianZ());
          beta=(manifoldY[k]->metric(rr,rr))/rrNorm;
          rrNorm=beta*rrNorm;
          if(rrNorm<rrNorm0){
            break;
          }
          dd=-rr+beta*dd;
        }
        //end of inner loop:: subproblem
        //move forward in the direction of eta with stepsize=1
          if(prodK>1){
             hessianF=as< arma::mat>(hessian(YList,k+1,eta));
          }else{
            hessianF=as< arma::mat>(hessian(YList[0],eta));
          }

        m_temp=manifoldY[k]->secondOrderApprox(objValue,eta,hessianF);

        manifoldY[k]->set_descD(eta);
       YList[k]=manifoldY[k]->retract(1,"trustRegion",true);
        if(prodK>1){
          objValue_temp=as< double>(obej(YList));
        }else{
          objValue_temp=as< double>(obej(YList[0]));
        }
        objDesc=objValue-objValue_temp;
        mDesc=objValue-m_temp;
        rho=objDesc/mDesc;
        if(rho<0.25){
           Delta[k]=0.25*Delta[k];
        }else if(rho>0.75 && bdryTouched){
          Delta[k]=min(2*Delta[k],DeltaMax);
        }
        if(rho>rhoMin){
          manifoldY[k]->acceptY();
          objValue=objValue_temp;
        }else{
          YList[k]=manifoldY[k]->get_Y();
        }
     }// iteration over product component
     objDesc=objValue_outer-objValue;
    //if(tol> DeltaL) flag=false;
   //  DeltaL=(*max_element(Delta.begin(),Delta.end()));
    if(rho>rhoMin*1.05 && objDesc<tol) flag=false;
    //if(rrNormMax<tol) flag=false;
  }// outer iteration

  return List::create(Named("optY")=YList,
              Named("optValue")=objValue,
              Named("NumIter")=iter);
END_RCPP
}//end of function
