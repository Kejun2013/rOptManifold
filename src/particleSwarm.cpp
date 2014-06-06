#include <omp.h>
#include "stiefel.h"
#include "grassmannQ.h"
#include <stdlib.h>
#include <time.h>


// YList is a (list of) of matrix(ces) of dimensional n1*p1, initial values
//f1 is objective function; f2 is gradient function
// iterMax: maximum number of iterations
// tol: tolerance of accuracy
//retract: retraction method

//// [[Rcpp::plugins(openmp)]]

RcppExport SEXP  particleSwarm(SEXP YList1, SEXP n1, SEXP p1, SEXP r1,
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
 /* double tol=as< double>(control["tol"]);
  double sigma=as< double>(control["sigma"]);
  double beta=as< double>(control["beta"]);
  double alpha=as< double>(control["alpha"]); */

  
  // Initialization of Data points
  IntegerVector n(n1),p(p1),r(r1);
  CharacterVector mtype(mtype1);
  int prodK=n.size();// size of product manifold
  vector< manifold*>  manifoldYG;
  int k;
  List YList(YList1), YList_temp(YList1);
  
  int dim=0;
  for(k=0;k<prodK;k++) dim+=n[k]+p[k];
  const int particle_num=dim; //dimension needs to be changed;
  vector< vector< manifold*> >  manifoldYP(particle_num), manifoldYB(particle_num);   
  

  
  int iter=0; // outer loop
  double objValue;
  vector<double> objValue_p(particle_num,0.0), objValue_b(particle_num,0.0);
    

 //initiating particle_num size of particles and velocities;
  int outter_num; 
  for(outter_num=0;outter_num<particle_num;outter_num++){
    for(k=0;k<prodK;k++){
      SEXP yTemp2=YList[k];
       NumericMatrix yTemp(yTemp2);
       std::string typeTemp=as< std::string>(mtype[k]);     
       if(typeTemp=="stiefel"){
          manifoldYP[outter_num].push_back(new stiefel(n[k],p[k],r[k],
                                    yTemp,retraction[k]));                   
          manifoldYB[outter_num].push_back(new stiefel(n[k],p[k],r[k],
                                    yTemp,retraction[k]));                        
          if(outter_num==0)  manifoldYG.push_back(new stiefel(n[k],p[k],r[k],
                                                     yTemp,retraction[k]));                          
       }
       else if(typeTemp=="grassmannQ"){
         manifoldYP[outter_num].push_back(new grassmannQ(n[k],p[k],r[k],
                                    yTemp,retraction[k]));
         manifoldYB[outter_num].push_back(new grassmannQ(n[k],p[k],r[k],
                                    yTemp,retraction[k]));     
          if(outter_num==0)  manifoldYG.push_back(new grassmannQ(n[k],p[k],r[k],
                                                     yTemp,retraction[k]));                                     
       }else if(typeTemp=="fixedRank"){
       
       }
       manifoldYP[outter_num][k]->set_particle();  
       *manifoldYB[outter_num][k]=*manifoldYP[outter_num][k];  
     //  if(outter_num==0)  *manifoldYG[k]=*manifoldYP[outter_num][k];
    }    
  } 
  
  
  
   
  //initialising the global-value postion
  for(k=0;k<prodK;k++) YList[k]=manifoldYG[k]->get_Y();
  if(prodK>1){
    objValue=as< double>(obej(YList));
  }else{
    objValue=as< double>(obej(YList[0]));
  }   
  for(outter_num=0;outter_num<particle_num;outter_num++){
    for(k=0;k<prodK;k++)  YList_temp[k]=manifoldYP[outter_num][k]->get_Y();
      if(prodK>1){
        objValue_b[outter_num]=as< double>(obej(YList_temp));
      }else{
        objValue_b[outter_num]=as< double>(obej(YList_temp[0]));
      }
      if(objValue_b[outter_num]<objValue) {
        for(k=0;k<prodK;k++) manifoldYG[k]->set_Y(manifoldYP[outter_num][k]->get_Y());
        objValue=objValue_b[outter_num];
      }  
  }

  //specific arguments of particleSwarm;
  double omega=1.0;  //can be added to arguments later;
  double phi1=2,phi2=2; //can be added to arguments later;
  
  arma::mat velocity;
  //R01 and R02 are uniform (0,1) distributed numbers;
  srand (time(NULL));
  double R01=0.5,R02=0.5;
  int thread_num;
  int subthread_num=0,subthread_num_1=0; //to test whether parallelism happens;
  
  //begin iteration  
  while(iter<iterMax){
    iter++;
    #pragma omp parallel shared(subthread_num,manifoldYG,thread_num,objValue,manifoldYB) \
       firstprivate(velocity,objValue_p,objValue_b,YList_temp) \
       private(R01,R02,k) 
    { 
      
      #pragma omp for schedule(static,8) lastprivate(subthread_num) ordered
        for(outter_num=0;outter_num<particle_num;outter_num++){
          if(omp_get_thread_num()==1) subthread_num_1++; //to test whether parallelism happens;
          thread_num=omp_get_num_threads();
          #pragma omp ordered
          {  
            for(k=0;k<prodK;k++){
              R01=((double) rand() / (RAND_MAX+1));
              R02=((double) rand() / (RAND_MAX+1));
              velocity=omega*manifoldYP[outter_num][k]->get_conjugateD();
              #pragma omp critical
              {
              velocity+=phi1*R01*(manifoldYB[outter_num][k]->get_Y()-
                                   manifoldYP[outter_num][k]->get_Y());
              velocity+=phi2*R02*(manifoldYG[k]->get_Y()-
                                    manifoldYP[outter_num][k]->get_Y());
              }  //out of critical                    
              manifoldYP[outter_num][k]->evalGradient(velocity,"steepest");
              velocity=manifoldYP[outter_num][k]->get_Gradient();
              YList_temp[k]=manifoldYP[outter_num][k]->genretract(1,velocity);
              velocity=manifoldYP[outter_num][k]->vectorTrans(1,velocity);
              manifoldYP[outter_num][k]->set_conjugateD(velocity); 
              manifoldYP[outter_num][k]->acceptY();
            }  
            if(prodK>1){
            objValue_p[outter_num]=as< double>(obej(YList_temp));
            }else{
            objValue_p[outter_num]=as< double>(obej(YList_temp[0]));
            }
            #pragma omp critical
            {
            if(objValue_p[outter_num]<objValue_b[outter_num]){ 
            objValue_b[outter_num]=objValue_p[outter_num];
            for(k=0;k<prodK;k++) *manifoldYB[outter_num][k]=*manifoldYP[outter_num][k];
            }  
            if(objValue_b[outter_num]<objValue) {
            for(k=0;k<prodK;k++) manifoldYG[k]->set_Y(manifoldYP[outter_num][k]->get_Y());
              objValue=objValue_b[outter_num];
            }
            }  //out of critical
          } //out of order          
        
         }// iteration over particles;
         #pragma omp single
         {
           subthread_num+=subthread_num_1;
         }
    } //out of parallel region;
  }// outer iteration

  for(k=0;k<prodK;k++) YList[k]=manifoldYG[k]->get_Y();  
  return List::create(Named("optY")=YList,
              Named("optValue")=objValue,
              Named("threads")=thread_num,
              Named("sub_threads")=subthread_num);
          //    Named("NumIter")=iter);
END_RCPP
}//end of function
