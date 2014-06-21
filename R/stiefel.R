stiefel<-function(n=10,p=2,obj=NULL,grad=NULL,hessian=NULL,retraction="QR"){
  if(length(n)!=length(p)) stop("N and P should have same length!")
  if(any(n<p)) stop("N should not be smaller than P!")
  if(length(n)==1){
    cat(paste0("Stiefel(",n,",",p,") Manifold Created!\n"))
  }else if(length(n)>1){
    message="Product Manifold Created: "
    for(i in 1:length(n)){
      message=paste0(message,"Stiefel(",n[i],",",p[i],")")
      if(i< length(n)) message=paste0(message, " * ") 
    }
    cat(message, "\n")
  }
  Y=lapply(1:length(n),function(iid){diag(1,n[iid],p[iid])})#initial value on manifold
  new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
      obj=obj,grad=grad,hessian=hessian,
      retraction=rep(retraction,length(n)),mtype=rep("stiefel",length(n))
  )
}


grassmanQ<-function(n=10,p=2,obj=NULL,grad=NULL,
                    hessian=NULL,retraction="Exp"){
  if(length(n)!=length(p)) stop("N and P should have same length!")
  if(any(n<p)) stop("N should not be smaller than P!")
  if(length(n)==1){
    cat(paste0("Grassman(",n,",",p,") Quotient Manifold Created!\n"))
  }else if(length(n)>1){
    message="Product Manifold Created: "
    for(i in 1:length(n)){
      message=paste0(message,"grassmanQ(",n[i],",",p[i],")")
      if(i< length(n)) message=paste0(message, " * ") 
    }
    cat(message, "\n")
  }
  Y=lapply(1:length(n),function(iid){diag(1,n[iid],p[iid])})#initial value on manifold
  new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
      obj=obj,grad=grad,hessian=hessian,
      retraction=rep(retraction,length(n)),
      mtype=rep("grassmanQ",length(n))
  )
}


fixedRank<-function(n=20,p=20,r=2,obj=NULL,grad=NULL,
                    hessian=NULL){
  if(length(n)!=length(p)) stop("N and P should have same length!")
  if(any(n<p)) stop("N should not be smaller than P!")
  if(length(n)==1){
    cat(paste0("Fixed Rank(",n,",",p,",",r,") Manifold Created!\n"))
  }else if(length(n)>1){
    message="Product Manifold Created: "
    for(i in 1:length(n)){
      message=paste0(message,"grassmanQ(",n[i],",",p[i],")")
      if(i< length(n)) message=paste0(message, " * ") 
    }
    cat(message, "\n")
  }
  Y=lapply(1:length(n),function(ii){
    temp=matrix(0,n[ii],p[ii])
    diag(temp)=c(rep(1,r),rep(0,min(n[ii],p[ii])-r))
    temp
  })#initial value on manifold
  new("manifold",Y=Y,n=n,p=p,r=r,
      obj=obj,grad=grad,hessian=hessian,
      retraction="QR",
      mtype=rep("fixedRank",length(n))
  )
}

#new added
sphere<-function(n=10,p=2,obj=NULL,grad=NULL,hessian=NULL,retraction="Exp"){
  if(length(n)!=length(p)) stop("N and P should have same length!")
#  if(any(n<p)) stop("N should not be smaller than P!")
  if(length(n)==1){
    cat(paste0("Sphere(",n,",",p,") Manifold Created!\n"))
  }else if(length(n)>1){
    message="Product Manifold Created: "
    for(i in 1:length(n)){
      message=paste0(message,"Sphere(",n[i],",",p[i],")")
      if(i< length(n)) message=paste0(message, " * ") 
    }
    cat(message, "\n")
  }
  Y=lapply(1:length(n),function(ii){
    temp=matrix(0,n[ii],p[ii])
    temp2=min(n[ii],p[ii])
    temp3=1/(sqrt(temp2))
    diag(temp)=c(rep(temp3,temp2))
    temp
  })#initial value on manifold
  new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
      obj=obj,grad=grad,hessian=hessian,
      retraction=rep(retraction,length(n)),mtype=rep("sphere",length(n))
  )
}




