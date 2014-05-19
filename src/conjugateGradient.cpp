//#include <omp.h>
#include "stiefel.h"
#include "grassmannQ.h"

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
  
  //besides the stepsize,we need one more eta as Fletcher-Reeves
  vector<double> eta;
  for (k=0;k<prodK;k++){  //captial K in prodK
    eta.push_back(0);
  }
  
  //
  
  //for conjugate gradient, we need conjugateD_temp;
  arma::mat conjugateD_temp;
  
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
  
  //need to start somewhere here 
  
  double gradNorm=1.0;
 // int iterInner=0; // to control the inner iterations;
  
  
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
        if(retraction[k]!=2) manifoldY[k]->evalGradient(gradF,"steepest");
        else manifoldY[k]->evalGradient(gradF,"a");
        if(eta[k]==0) manifoldY[k]->set_conjugateD(-(manifoldY[k]->get_Gradient()));
        
        
      //  gradNorm=arma::dot(manifoldY[k]->get_Gradient(),manifoldY[k]->get_Gradient());
        gradNorm=manifoldY[k]->metric(manifoldY[k]->get_Gradient(),manifoldY[k]->get_Gradient());
        stepsize=alpha/beta;
        eDescent=sigma/beta*manifoldY[k]->metric(manifoldY[k]->get_Gradient(),-(manifoldY[k]->get_conjugateD()));
        do{//c
          stepsize=stepsize*beta;
          eDescent=eDescent*beta;
          YList[k]=manifoldY[k]->genretract(stepsize,manifoldY[k]->get_conjugateD());
          if(prodK>1){
            objValue_temp=as< double>(obej(YList));
          }else{
            objValue_temp=as< double>(obej(YList[0]));
          }
        }while((objValue-objValue_temp)<eDescent);
        //step size iteration
        //if a stepsize is accepted, update current location
        conjugateD_temp= manifoldY[k]->vectorTrans(stepsize,manifoldY[k]->get_conjugateD());   
        manifoldY[k]->acceptY();   
        //update gradient and conjudate gradient
        if(retraction[k]!=2) manifoldY[k]->evalGradient(gradF,"steepest");
        else manifoldY[k]->evalGradient(gradF,"a");
        
        eta[k]=manifoldY[k]->metric(manifoldY[k]->get_Gradient(),manifoldY[k]->get_Gradient())/(gradNorm);
        conjugateD_temp*=eta[k];
        conjugateD_temp+=-(manifoldY[k]->get_Gradient());
        
        manifoldY[k]->set_conjugateD(conjugateD_temp);        
        objDesc=objValue-objValue_temp;
        objValue=objValue_temp;
    }// iteration over product component
    if(tol>objDesc) { 
      double tolInner=0.0;  // in CG, sometimes conjugate=0,while grad!=0
      for(k=0;k<prodK;k++) tolInner+=manifoldY[k]->metric(manifoldY[k]->get_Gradient(),manifoldY[k]->get_Gradient());
      if(abs(objValue)*tol>sqrt(tolInner)) flag=false;
     // else eta[k]=0;
    }
  }// outer iteration

  return List::create(Named("optY")=YList,
              Named("optValue")=objValue,
              Named("NumIter")=iter);
END_RCPP
}//end of function




