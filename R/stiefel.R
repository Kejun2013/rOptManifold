setClassUnion("NullFunction",c("NULL","function"))

setClass("stiefel",
        representation(n="numeric",p="numeric",obj="NullFunction",
                        grad="NullFunction",iterMax="numeric",
                       tol="numeric",
                       retraction="character"),
         prototype=prototype(n=10,p=2,obj=NULL,grad=NULL,iterMax=100,
                             retraction="QR"))

stiefel<-function(n=10,p=2,obj=NULL,grad=NULL,retraction="QR"){
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
  new("stiefel",n=n,p=p,obj=obj,grad=grad,retraction=retraction,tol=1e-6)
}

setGeneric("setObjective<-",function(object,value){standardGeneric("setObjective<-")})
setGeneric("setGradient<-",function(object,value){standardGeneric("setGradient<-")})
setGeneric("steepestDescent",function(object){standardGeneric("steepestDescent")})


setReplaceMethod(
  f="setObjective",
  signature="stiefel",
  definition=function(object,value){
    object@obj<-value
    return(object)
  }
)

setReplaceMethod(
  f="setGradient",
  signature="stiefel",
  definition=function(object,value){
    object@grad<-value
    return(object)
  }
)

setMethod(
  f="[",
  signature="stiefel",
  definition=function(x,i,j,drop){
    if(i=="obj" || i=="objective"){return(x@obj)}
    if(i=="grad" || i=="gradient"){return(x@grad)}
    if(i=="iterMax"){return(x@iterMax)}
    if(i=="retract" || i=="retraction"){return(x@retraction)}
    if(i=="tol"){return(x@tol)}
  }
)

setReplaceMethod(
  f="[",
  signature="stiefel",
  definition=function(x,i,j,value){
    if(i=="obj" || i=="objective"){x@obj<-value}
    if(i=="grad" || i=="gradient"){x@grad<-value}
    if(i=="iterMax"){x@iterMax<-value}
    if(i=="retract" || i=="retraction"){x@retraction<-value}
    if(i=="tol"){x@tol<-value}
    return(x)
  }
)

setMethod("show","stiefel",
          function(object){
            if(length(object@n)==1){
              cat(paste0("Stiefel(",object@n,",",object@p,") Manifold"))
            }else if(length(object@n)>1){
              message="Product Manifold: "
              for(i in 1:length(object@n)){
                message=paste0(message,"Stiefel(",object@n[i],",",object@p[i],")")
                if(i< length(object@n)) message=paste0(message, " * ") 
              }
              cat(message)
            }
            if(is.null(object@obj)){
              cat("\n    Objective Function: Not Available")
            }else{
              cat("\n    Objective Function: Available")
            }
            if(is.null(object@grad)){
              cat("\n    Gradient Function: Not Available")
            }else{
              cat("\n    Gradient Function: Available")
            }
            cat("\n    Retraction Method:", toupper(object@retraction))
          })

setMethod("*",signature(e1="stiefel",e2="stiefel"),
          function(e1,e2){
            n=c(e1@n,e2@n)
            p=c(e1@p,e2@p)
            stiefel(n,p)
          })

setMethod("steepestDescent","stiefel",
          definition=function(object){
            retractMethod=switch(tolower(object@retraction),exp=0,qr=1,cayley=2)
            .Call("steepestDescent_stiefel",
                  object@n,object@p,object@obj,object@grad,
                  function(X){as.matrix(expm(X))},
                  object@iterMax,object@tol,retractMethod,
                  PACKAGE = "rOptManifold" )
          })

