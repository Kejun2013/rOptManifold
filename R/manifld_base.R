####This fle contains basic definition of Manfiold Class and related operation


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
                             control=list(iterMax=1000,###common
                                          iterSubMax=50,
                                          tol=0.0001,
                                          conjMethod="PR",##conj and steepest
                                          alpha=5,beta=0.8,sigma=0.6,
                                          theta=1,kappa=0.01,rhoMin=0.1,##trustRegion
                                          Delta0=0.5,DeltaMax=5,
                                          phi1=2,phi2=2,omega=1,##para particle swarm
                                          particleNum=50,threadNum=4)
         )
)


setGeneric("setGradient<-",function(object,value){standardGeneric("setGradient<-")})



setReplaceMethod(
  f="setGradient",
  signature="manifold",
  definition=function(object,value){
    object@grad<-value
    return(object)
  }
)

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
