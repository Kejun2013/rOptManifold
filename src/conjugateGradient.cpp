//#include <omp.h>
#include "stiefel.h"
#include "grassmannQ.h"
#include "fixRank.h"
#include "fixRankPSD.h"
#include "fixRankSym.h"
#include "spectahedron.h"
#include "elliptope.h"
#include "sphere.h"
#include "oblique.h"
#include "specialLinear.h"
#include "projective.h"
//' Conjugate gradient method
//' 
//' @param YList a (list of) of matrix(ces) of dimensional n1*p1, initial values
//' @param f1 objective function
//' @param f2 gradient function
//' @param iterMax in control, maximum number of iterations
//' @param tol in control, tolerance of accuracy
//' @param alpha in control, the initial stepsize 
//' @param beta in control, beta 
//' @param sigma in control, Armijo stepsize controllers
//' @param retraction1 retraction method represented by int
//' @return an R list of argmin, mininum values and number of iterations 


//The conjugate direction is stored in descD;

RcppExport SEXP  conjugateGradient(SEXP YList1, SEXP n1, SEXP p1, SEXP r1,
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
  //0: Fletcher-Reeves 1: Polak-Ribiere

   string conjMethod1=as< string>(control["conjMethod"]);
   int conjMethod;
   if(conjMethod1=="PR") conjMethod=1;
   else conjMethod=0;
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
     }else if(typeTemp=="specialLinear"){
       manifoldY.push_back(new specialLinear(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }else if(typeTemp=="projective"){
       manifoldY.push_back(new projective(n[k],p[k],r[k],
                                  yTemp,retraction[k]));
     }       

  }
    
  //define other varibles
  int iter=0,iterInner=0; // outer loop, inner loop control
  bool flag=true,first=true;
  double objValue,objValue_outer,objValue_temp,objDesc,stepsize;
  double gradNorm_temp=0,ang_check=0;
  //besides the stepsize,we need one more eta as Fletcher-Reeves
  vector< double> gradNorm(prodK,1.0),eDescent(prodK,1.0);
  
  
  arma::mat gradF,xi_bar;  //gradient in ambient space
  
  if(prodK>1){
    objValue=as< double>(obej(YList));
  }else{
    objValue=as< double>(obej(YList[0]));
  }
  
  //need to start somewhere here 
  
  double eta=0.0;
 // int iterInner=0; // to control the inner iterations;
  
  
  //initializing conjugate direction
  for(k=0;k<prodK;k++){
      if(prodK>1){
        gradF=as< arma::mat>(grad(YList,k+1));
      }else{
        gradF=as< arma::mat>(grad(YList[0]));
      }
        //gradient on the stiefel manifold, set descD to -xi
      manifoldY[k]->evalGradient(gradF,"steepest");
      gradNorm[k]=manifoldY[k]->gradMetric(); //<gradF,gradF>_Y
      eDescent[k]=gradNorm[k];
  }
  
  while(iter<iterMax && flag){
    //gradient of objective funtion
    iter++;
    objValue_outer=objValue;
    for(k=0;k<prodK;k++){
        stepsize=alpha/beta;
        //Rcpp::Rcout<<"eD="<<eDescent[k]<<endl;
        eDescent[k]=sigma/beta*(eDescent[k]); //<descD,grad F>_Y
        eDescent[k]=abs(eDescent[k]);
        first=true;
        iterInner=0;
        //manifoldY[k]->set_descD(manifoldY[k]->get_conjugateD());
        
        do{//choose appropriate stepsize by Armijo rule
          iterInner++;
          stepsize=stepsize*beta;
          eDescent[k]=eDescent[k]*beta;
          YList[k]=manifoldY[k]->retract(stepsize,"conjugateD",first);
          if(prodK>1){
            objValue_temp=as< double>(obej(YList));
          }else{
            objValue_temp=as< double>(obej(YList[0]));
          }
          first=false;
        }while((objValue-objValue_temp)<eDescent[k] && iterInner<1000);
        //end of step size iteration
        
        
        //vector trans of descD to Yt
          manifoldY[k]->vectorTrans();
         xi_bar=manifoldY[k]->get_Gradient(); //current gradient
         
         //move toward new position and compute next descD
          manifoldY[k]->acceptY();
          objValue=objValue_temp;
          
          //compute new gradient at new position
          if(prodK>1){
            gradF=as< arma::mat>(grad(YList,k+1));
          }else{
            gradF=as< arma::mat>(grad(YList[0]));
          }
          manifoldY[k]->evalGradient(gradF,"conjugate");
          
          //decide eta for descD_new=-xi_Yt+eta*descD_transported
          if(conjMethod==0){//Fletcher-Reeves
            eta=(manifoldY[k]->gradMetric())/gradNorm[k];
            gradNorm[k]=eta*gradNorm[k];//update gradient norm at current position
          }else{//Polak-Ribiere, conjMethod==1
            gradNorm_temp=manifoldY[k]->gradMetric();
            eta=(gradNorm_temp-
                    manifoldY[k]->metric(xi_bar,
                                manifoldY[k]->get_Gradient()))/gradNorm[k];
            eta=std::max(0.0,eta);
            gradNorm[k]=gradNorm_temp;
          }
          //update descD according to descD_new=-xi_Yt+eta*descD_transported
          //done after transport descD and compute new gradient
          manifoldY[k]->update_conjugateD(eta);
          
          //check the angel of xi and descD, if >-(0.05) set descD=-xi
           //this makes sure descD is a descent direction
           ang_check=manifoldY[k]->grad_descD_Metric();
           eDescent[k]=ang_check;
           ang_check=ang_check/sqrt(gradNorm[k]);
           ang_check=ang_check/sqrt(manifoldY[k]->descDMetric());
           
           if(ang_check>(-0.05)){
                Rcpp::Rcout<<"Iteration "<<iter<<":"<<endl;
                Rcpp::Rcout<<"Conjugate Direction Reset To Steepest Descent!"<<endl;
                manifoldY[k]->retrieve_steepest();
                eDescent[k]=gradNorm[k];
           }
          //after moving all directions from Y to Yt, set Y=Yt

    }// iteration over product component
    objDesc=objValue_outer-objValue;
    //Rcpp::Rcout<<"desc="<<objDesc<<endl;
    if(tol>objDesc) { 
      flag=false;
//      double tolInner=0.0;  // in CG, sometimes conjugate=0,while grad!=0
//      for(k=0;k<prodK;k++) tolInner+=gradNorm[k];
//      if(abs(objValue)*tol>sqrt(tolInner)) flag=false;
//     // else eta[k]=0;
    }
  }// outer iteration

  return List::create(Named("optY")=YList,
              Named("optValue")=objValue,
              Named("NumIter")=iter);
END_RCPP
}//end of function




