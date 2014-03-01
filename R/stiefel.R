setClassUnion("NullFunction",c("NULL","function"))

setClass("stiefel",
        representation(n="numeric",p="numeric",obj="NullFunction",
                        grad="NullFunction",iterMax="numeric",retraction="character"),
         prototype=prototype(n=10,p=2,obj=NULL,grad=NULL,iterMax=100,
                             retraction="QR"))

stiefel<-function(n=10,p=2,obj=NULL,grad=NULL){
  new("stiefel",n=n,p=2,obj=obj,grad=grad)
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
    return(x)
  }
)


setMethod("steepestDescent","stiefel",
          definition=function(object){
            retractMethod=switch(tolower(object@retraction),exp=0,qr=1,cayley=2)
            .Call("rcpparma_hello_world",
                  object@n,object@p,object@obj,object@grad,
                  function(X){as.matrix(expm(X))},
                  object@iterMax,retractMethod,
                  PACKAGE = "rOptManifold" )
          })

