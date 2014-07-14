########Optimization Methods in this packages
########Including: (1) Steepest Descent (2) Trust Region (3) Conjugate Gradient (4) Particle Swarm


###This function check user specified retraction method
###If user mis-specified, set to default
SetRetraction=function(object){
  retractMethod=rep(0,length(object@n))
  for(i in 1:length(object@n)){
    retr=tolower(object@retraction[i])
    
    if(object@mtype[i]=="stiefel"){
        if(!(retr %in% c("qr","cayley"))){
          cat(paste0("Component ",i,": ",object@mtype[i],
                        " retraction should be one of 'QR', 
                        'Cayley'. Set to 'QR' instead"))
          object@retraction[i]="QR" 
        }
    }else if(object@mtype[i]=="sphere"){
       if(!(retr %in% c("exp","norm"))){
         cat(paste0("Component ",i,": ",object@mtype[i],
                        " retraction should be one of 'Exp', 'Norm'. 
                        Set to 'Norm' instead"))
         object@retraction[i]="Norm" 
       }
    }else if(object@mtype[i]=="grassmannQ"){
      if(!(retr %in% c("exp","qr"))){
        cat(paste0("Component ",i,": ",object@mtype[i],
                       " retraction should be one of 'QR', 'Exp'. 
                       Set to 'QR' instead"))
        object@retraction[i]="QR" 
      }
    }else if(object@mtype[i]=="grassmannSub"){
      if(!(retr %in% c("qr","cayley"))){
        cat(paste0("Component ",i,": ",object@mtype[i],
                       " retraction should be one of 'QR', 'Cayley'.
                       Set to 'QR' instead"))
        object@retraction[i]="QR" 
      }
    }else if(object@mtype[i] %in% c("fixedRank","fixedRankSym")){
      if(!(retr %in% c("proj"))){
        cat(paste0("Component ",i,": ",object@mtype[i],
                       " retraction should be  'Proj'. 
                       Set to 'Proj' instead"))
        object@retraction[i]="Proj" 
      }
    }else if(object@mtype[i] %in% c("fixedRankPSD","spectahedron","elliptope")){
      if(!(retr %in% c("exp","proj"))){
        cat(paste0("Component ",i,": ",object@mtype[i],
                   " retraction should be one of 'Exp', 'Proj'. 
                   Set to 'Proj' instead"))
        object@retraction[i]="Proj" 
      }
    }
    
    else if(object@mtype[i]=="oblique"){
      if(!(retr %in% c("norm"))){
        cat(paste0("Component ",i,": ",object@mtype[i],
                   " retraction should be  'Norm'. 
                   Set to 'Norm' instead"))
        object@retraction[i]="Norm" 
      }
    }else{
      stop(paste("Manifold Type Not Found:", object@mtype[i] ))
    }
    retractMethod[i]=switch(tolower(object@retraction[i]),exp=0,
                                                     qr=1,cayley=2,
                                                     proj=3,norm=1) 
  }
  
  retractMethod
}


store_obj=function(){}
store_grad=function(){}
num_n=1
num_p=1
obj2grad=function(Y,k=0){
  numericalGradient(store_obj,Y,num_n,num_p,k)
}
obj2hess=function(Y,Z,k=0){
  numericalHessian(store_obj,Y,Z,num_n,num_p,k)
}
grad2hess=function(Y,Z,k=0){
  numericalHessian2(store_grad,Y,Z,k)
}




setGeneric("steepestDescent",function(object){standardGeneric("steepestDescent")})
setMethod("steepestDescent","manifold",
          definition=function(object){
            retractMethod=SetRetraction(object)
            optUse_grad=object@grad
            if(is.null(optUse_grad)){
              assignInNamespace("store_obj",object@obj,"rOptManifold")
              assignInNamespace("num_n",object@n,"rOptManifold")
              assignInNamespace("num_p",object@p,"rOptManifold")
              optUse_grad=obj2grad
            }
            
            .Call("steepestDescent",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,optUse_grad,
                  object@control,
                  PACKAGE = "rOptManifold" )
          })



setGeneric("trustRegion",function(object){standardGeneric("trustRegion")})
setMethod("trustRegion","manifold",
          definition=function(object){
            retractMethod=SetRetraction(object)
            optUse_grad=object@grad
            assignInNamespace("store_obj",object@obj,"rOptManifold")
            assignInNamespace("store_grad",object@grad,"rOptManifold")
            assignInNamespace("num_n",object@n,"rOptManifold")
            assignInNamespace("num_p",object@p,"rOptManifold")
            if(is.null(optUse_grad)){
              optUse_grad=obj2grad
            }
            optUse_hess=object@hessian
            if(is.null(optUse_hess)){
              if(is.null(object@grad)){
                optUse_hess=obj2hess
              }else{
                optUse_hess=grad2hess
              }
            }
            .Call("trustRegion",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,optUse_grad,optUse_hess,
                  object@control,
                  PACKAGE = "rOptManifold" )
          })


## Kejun just does not know how to write this one
setGeneric("conjugateGradient",function(object){standardGeneric("conjugateGradient")})
setMethod("conjugateGradient","manifold",
          definition=function(object){
            retractMethod=SetRetraction(object)
            optUse_grad=object@grad
            if(is.null(optUse_grad)){
              assignInNamespace("store_obj",object@obj,"rOptManifold")
              assignInNamespace("num_n",object@n,"rOptManifold")
              assignInNamespace("num_p",object@p,"rOptManifold")
              optUse_grad=obj2grad
            }
            .Call("conjugateGradient",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,optUse_grad,
                  object@control,
                  PACKAGE = "rOptManifold")
          })

## Kejun just does not know how to write this one
setGeneric("particleSwarm",function(object){standardGeneric("particleSwarm")})
setMethod("particleSwarm","manifold",
          definition=function(object){
            retractMethod=SetRetraction(object)
            .Call("particleSwarm",
                  object@Y,
                  object@n,object@p,object@r,
                  object@mtype,retractMethod, 
                  object@obj,
                  object@control,
                  PACKAGE = "rOptManifold")
          })
