setClassUnion("NullFunction",c("NULL","function"))


## n-by-p matrix of rank r (r is required only for some type of manifolds)
setClass("manifold",
         representation(Y="list",
                        n="numeric",p="numeric",r="numeric",
                        mtype="character",
                        obj="NullFunction",
                        grad="NullFunction",
                        hessian="NullFunction",
                        retraction="character",
                        control="list"),
         prototype=prototype(Y=list(diag(1,10,2)),
                            n=10,p=2,r=0,
                             mtype="stiefel",
                             obj=NULL,grad=NULL,hessian=NULL,
                             retraction="QR",
                             control=list(iterMax=1000,tol=0.0001,
                                          alpha=5,beta=0.8,sigma=0.6,
                                          theta=1,kappa=0.01,rhoMin=0.1,
                                          Delta0=0.5,
                                          DeltaMax=5)
                             )
         )


# setGeneric("setGradient<-",function(object,value){standardGeneric("setGradient<-")})


# 
# setReplaceMethod(
#   f="setGradient",
#   signature="manifold",
#   definition=function(object,value){
#     object@grad<-value
#     return(object)
#   }
# )

setMethod(
  f="[",
  signature="manifold",
  definition=function(x,i,j,drop){
    if(i=="obj" || i=="objective"){return(x@obj)}
    if(i=="grad" || i=="gradient"){return(x@grad)}
    if(i=="hessian"){return(x@hessian)}
    if(i=="retract" || i=="retraction"){return(x@retraction[j])}
    if(i=="control" ){return(x@control)}
  }
)

setReplaceMethod(
  f="[",
  signature="manifold",
  definition=function(x,i,j,value){
    i=tolower(i)
    if(i=="obj" || i=="o"){x@obj<-value}
    if(i=="grad" || i=="g"){x@grad<-value}
    if(i=="hessian" || i=="h"){x@hessian<-value}
    if(i=="retraction" || i=="r"){x@retraction<-value} # NOt retractioin
    if(i=="control" || i=="c"){x@control[[j]]<-value}
    return(x)
  }
)

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



# 
# setMethod("show","manifold",
#           function(object){
#             if(length(object@n)==1){
#               cat(paste0("Stiefel(",object@n,",",object@p,") Manifold"))
#             }else if(length(object@n)>1){
#               message="Product Manifold: "
#               for(i in 1:length(object@n)){
#                 message=paste0(message,"Stiefel(",object@n[i],",",object@p[i],")")
#                 if(i< length(object@n)) message=paste0(message, " * ") 
#               }
#               cat(message)
#             }
#             if(is.null(object@obj)){
#               cat("\n    Objective Function: Not Available")
#             }else{
#               cat("\n    Objective Function: Available")
#             }
#             if(is.null(object@grad)){
#               cat("\n    Gradient Function: Not Available")
#             }else{
#               cat("\n    Gradient Function: Available")
#             }
#             cat("\n    Retraction Method:", toupper(object@retraction))
#           })

setMethod("*",signature(e1="manifold",e2="manifold"),
          function(e1,e2){
            Y=c(e1@Y,e2@Y)
            n=c(e1@n,e2@n)
            p=c(e1@p,e2@p)
            r=c(e1@r,e2@r)
            mtype=c(e1@mtype,e2@mtype)
            retraction=c(e1@retraction,e2@retraction)
            new("manifold",Y=Y,n=n,p=p,r=r,
                obj=NULL,grad=NULL,hessian=NULL,
                retraction=retraction,mtype=mtype)
          })


setGeneric("steepestDescent",function(object){standardGeneric("steepestDescent")})
setMethod("steepestDescent","manifold",
          definition=function(object){
            retractMethod=rep(0,length(object@n))
            for(i in 1:length(object@n)){
              retractMethod[i]=switch(tolower(object@retraction[i]),exp=0,qr=1,cayley=2) 
            }
            .Call("steepestDescent",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,object@grad,
                  object@control,
                  PACKAGE = "rOptManifold" )
          })



setGeneric("trustRegion",function(object){standardGeneric("trustRegion")})
setMethod("trustRegion","manifold",
          definition=function(object){
            retractMethod=rep(0,length(object@n))
            for(i in 1:length(object@n)){
              retractMethod[i]=switch(tolower(object@retraction[i]),exp=0,qr=1,cayley=2) 
            }
            .Call("trustRegion",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,object@grad,object@hessian,
                  object@control,
                  PACKAGE = "rOptManifold" )
          })


## Kejun just does not know how to write this one
setGeneric("conjugateGradient",function(object){standardGeneric("conjugateGradient")})
setMethod("conjugateGradient","manifold",
          definition=function(object){
            retractMethod=rep(0,length(object@n))
            for(i in 1:length(object@n)){
              retractMethod[i]=switch(tolower(object@retraction[i]),"exp"=0,"qr"=1,"cayley"=2) 
            }
            .Call("conjugateGradient",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,object@grad,
                  object@control,
                  PACKAGE = "rOptManifold" )
          })

## Kejun just does not know how to write this one
setGeneric("particleSwarm",function(object){standardGeneric("particleSwarm")})
setMethod("particleSwarm","manifold",
          definition=function(object){
            retractMethod=rep(0,length(object@n))
            for(i in 1:length(object@n)){
              retractMethod[i]=switch(tolower(object@retraction[i]),"exp"=0,"qr"=1,"cayley"=2) 
            }
            .Call("particleSwarm",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,object@grad,
                  object@control,
                  PACKAGE = "rOptManifold" )
          })

