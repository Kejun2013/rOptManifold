######functions for creating manifolds#####

####Basis Consistency Check for manifolds#######
BasicCheck<-function(n,p,r,retraction){
  if(length(n)!=length(p)) stop("N and P should have same length!")
  if(!is.null(r)){
    if(length(n)!=length(r) || length(p)!=length(r)) 
      stop("N, P and R should have same length!")
    if(any(r>p) || any(r>n)) 
      stop("R should not be greater than N or P!")
  }
  if(length(retraction)==1){
    retraction=rep(retraction,length(n))
  }else if(length(retraction)!=length(n)){
    stop("Retraction should have length of one, or same length as N!")
  }
  retraction
}


outputMessage<-function(n,p,r,mtype){
  if(is.null(r)){
    if(length(n)==1){
      cat(paste0(mtype,"(n=", n,",p=", p,") Manifold Created!\n"))
    }else if(length(n)>1){
      message="Product Manifold Created: "
      for(i in 1:length(n)){
        message=paste0(message,mtype,"(",n[i],",",p[i],")")
        if(i< length(n)) message=paste0(message, " * ") 
      }
      cat(message, "\n")
    }
  }else{
    if(length(n)==1){
      cat(paste0(mtype,"(n=", n,",p=", p,",r=",r,") Manifold Created!\n"))
    }else if(length(n)>1){
      message="Product Manifold Created: "
      for(i in 1:length(n)){
        message=paste0(message,mtype,"(",n[i],",",p[i],",", r[i] ,")")
        if(i< length(n)) message=paste0(message, " * ") 
      }
      cat(message, "\n")
    }
  }
}


#new added
sphere<-function(n,p,retraction="Norm"){
  retraction=BasicCheck(n,p,NULL,retraction)
  Y=lapply(1:length(n),function(ii){
    temp=matrix(0,n[ii],p[ii])
    diag(temp)=c(1,rep(0,min(n[ii],p[ii])-1))
    temp
  })#initial value on manifold
  mtype=rep("sphere",length(n))
  
  object=new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
      retraction=retraction,mtype=mtype)
  outputMessage(n,p,NULL,mtype)
  object
}

stiefel<-function(n,p,retraction="QR"){
  retraction=BasicCheck(n,p,NULL,retraction)
  if(any(n<p)) stop("N should not be smaller than P!")
  
  Y=lapply(1:length(n),function(iid){diag(1,n[iid],p[iid])})#initial value on manifold
  mtype=rep("stiefel",length(n))
  
  object=new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
      retraction=retraction,mtype=mtype)
  outputMessage(n,p,NULL,mtype)
  object
}


grassmanQ<-function(n,p,retraction="Exp"){
  retraction=BasicCheck(n,p,NULL,retraction)
  if(any(n<p)) stop("N should not be smaller than P!")

  Y=lapply(1:length(n),function(iid){diag(1,n[iid],p[iid])})#initial value on manifold
  mtype=rep("grassmanQ",length(n))
  
  object=new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
      retraction=retraction,mtype=mtype)
  
  outputMessage(n,p,NULL,mtype)
  object
}


fixedRank<-function(n,p,r, retraction="Proj"){
  retraction=BasicCheck(n,p,r,retraction)
  
  Y=lapply(1:length(n),function(ii){
    temp=matrix(0,n[ii],p[ii])
    diag(temp)=c(rep(1,r[ii]),rep(0,min(n[ii],p[ii])-r[ii]))
    temp
  })#initial value on manifold
  
  mtype=rep("fixedRank",length(n))
  object=new("manifold",Y=Y,n=n,p=p,r=r,
      retraction=retraction,mtype=mtype)
  outputMessage(n,p,r,mtype)
  object
}



fixedRankSym<-function(n,r,retraction="Proj"){
  retraction=BasicCheck(n,n,r,retraction)
  Y=lapply(1:length(n),function(ii){
    temp=matrix(0,n[ii],n[ii])
    diag(temp)=c(rep(1,r[ii]),rep(0,n[ii]-r[ii]))
    temp
  })#initial value on manifold
  mtype=rep("fixedRankSym",length(n))
  object=new("manifold",Y=Y,n=n,p=n,r=r,
      retraction=retraction,mtype=mtype)
  outputMessage(n,n,r,mtype)
  object
}



fixedRankPSD<-function(n,r,retraction="Proj"){
  retraction=BasicCheck(n,n,r,retraction)
  Y=lapply(1:length(n),function(ii){
    diag(1,n[ii],r[ii])
  })#initial value on manifold
  mtype=rep("fixedRankPSD",length(n))
  object=new("manifold",Y=Y,n=n,p=r,r=r,
      retraction=retraction,mtype=mtype)
  outputMessage(n,n,r,mtype)
  object
}


spectahedron<-function(n,r,retraction="Proj"){
  retraction=BasicCheck(n,n,r,retraction)
  Y=lapply(1:length(n),function(ii){
    temp=diag(1/r[ii],n[ii],r[ii])
  })#initial value on manifold
  mtype=rep("spectahedron",length(n))
  object=new("manifold",Y=Y,n=n,p=r,r=r,
      retraction=retraction,mtype=mtype)
  outputMessage(n,n,r,mtype)
  object
}



elliptope<-function(n,r,retraction="Proj"){
  retraction=BasicCheck(n,n,r,retraction)
  Y=lapply(1:length(n),function(ii){
    temp=diag(1,n[ii],r[ii])
    temp[(r[ii]+1):(n[ii]),1]=1
    temp
  })#initial value on manifold
  mtype=rep("elliptope",length(n))
  object=new("manifold",Y=Y,n=n,p=r,r=r,
      retraction=retraction,mtype=mtype)
  outputMessage(n,n,r,mtype)
  object
}

oblique<-function(n,p,retraction="Norm"){
  retraction=BasicCheck(n,p,NULL,retraction)
  Y=lapply(1:length(n),function(ii){
    diag(1,n[ii],p[ii])
  })#initial value on manifold
  mtype=rep("oblique",length(n))
  
  object=new("manifold",Y=Y,n=n,p=p,r=rep(0,length(n)),
             retraction=retraction,mtype=mtype)
  outputMessage(n,p,NULL,mtype)
  object
}
