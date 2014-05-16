//#include <omp.h>
#include "stiefel.h"
#include "grassmannQ.h"

//' Conjugate gradient method
//' 
//'  @param YList is a (list of) of matrix(ces) of dimensional n1*p1, initial values
//' @param f1 is objective function; f2 is gradient function
//' @param in control, iterMax: maximum number of iterations
//' @param in control, tol: tolerance of accuracy
//' @param in control, alpha: the initial stepsize 
//' @param in control, beta and sigma: Armijo stepsize controllers
//' @param retraction1: retraction method represented by int
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
  arma::mat conjugateD_temp,conjugate_pre;
  
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
        
        //conjugate_pre=(eta[k]==0)? -(manifoldY[k]->get_Gradient()): manifold[k]->get_conjugate();
       // conjugate_pre=manifold[k]->get_conjugate();
       if(eta[k]==0) {
         manifoldY[k]->evalGradient(gradF,"steepest");
         manifoldY[k]->set_conjugateD(-(manifoldY[k]->get_Gradient()));
       }
        
        conjugate_pre=manifoldY[k]->get_conjugateD();
        gradNorm=arma::dot(manifoldY[k]->get_Gradient(),manifoldY[k]->get_Gradient());
        stepsize=alpha/beta;
        
        eDescent=sigma/beta*(arma::dot(manifoldY[k]->get_Gradient(),-conjugate_pre));
        do{//c
          stepsize=stepsize*beta;
          eDescent=eDescent*beta;
          YList[k]=manifoldY[k]->genretract(stepsize,conjugate_pre);
          if(prodK>1){
            objValue_temp=as< double>(obej(YList));
          }else{
            objValue_temp=as< double>(obej(YList[0]));
          }
        }while((objValue-objValue_temp)<eDescent);
        //step size iteration
        //if a stepsize is accepted, update current location
        
       // eta[k]=arma::dot(YList[k],YList[k])/(arma::dot(manifold[k]->get_Y(),manifold[k]->get_Y()); 
        conjugateD_temp= manifoldY[k]->vectorTrans(stepsize,conjugate_pre);                  
        manifoldY[k]->acceptY();   
        manifoldY[k]->evalGradient(gradF,"steepest");   
        eta[k]=arma::dot(manifoldY[k]->get_Gradient(),manifoldY[k]->get_Gradient())/(gradNorm);
        
        conjugateD_temp*=eta[k];
        conjugateD_temp+=-(manifoldY[k]->get_Gradient());
        
        manifoldY[k]->set_conjugateD(conjugateD_temp);
        
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




