#include "manifold.h"

double manifold::get_eDescent(){
  return eDescent;
}

void manifold::acceptY(){
  Y=Yt;
}


manifold::manifold(int n1, int p1, int r1, 
              Rcpp::NumericMatrix Y1,int retraction1){
  n=n1;
  p=p1;
  r=r1;
  Y=arma::mat(Y1.begin(),n,p,false);
  retraction=retraction1;
}

