n=40##even number
set.seed(50)
B=C=D=matrix(0,n/2,n/2)

##The graph is devided into two group
##First group is 1:(n/2), the second groupd is (n/2+1):n
##Within each group, with is 0.2 prob that an edge exsits
##Between group, the prob is 0.6
fHalf=1:(n/2);sHalf=(n/2+1):n
elSel=upper.tri(B)
B[elSel]=(runif(sum(elSel))<0.03);
C[elSel]=(runif(sum(elSel))<0.3)
D[elSel]=(runif(sum(elSel))<0.05)

##A is adjacency matrix and L is Laplacian matrix
A=matrix(0,n,n)
A[fHalf,fHalf]=B;
A[fHalf,sHalf]=C;
A[sHalf,sHalf]=D;
A=A+t(A)
w=colSums(A)
L=diag(w)-A

par(mar=rep(0,4))
plot(c(1,5),c(1,(n/4+1)),type="n",bty="n",xaxt="n",yaxt="n")
px=c(rep(1.6,n/4),rep(2,n/4),rep(4,n/4),rep(4.6,n/4))
py=c(1:(n/4)+0.5,1:(n/4),1:(n/4)+0.5,1:(n/4))

for(i in 1:n){
  for(j in 1:n){
    if(A[i,j]==1){
      segments(px[i],py[i],px[j],py[j])
    }
  }
}

for(i in 1:n){
  points(px[i],py[i],pch=16,col="grey",cex=3)
}


rr=6
problem=elliptope(n=n,r=rr)
problem["obj"]=function(Y){
  objValue=sum(L*(Y%*%t(Y)))
  -objValue
}
#set gradient function
problem["grad"]=function(Y){
  -2*L%*%Y
}
#set hessian function
problem["hessian"]=function(Y,Z){
  -2*L%*%Z
}
problem["retraction"]="proj"
problem["control","tol"]=0.1
problem["control","Delta0"]=1
problem["control","DeltaMax"]=5
problem["control","rhoMin"]=0.1
problem["control","iterMax"]=1000
problem["control","iterSubMax"]=2000
problem["control","alpha"]=100

#res=trustRegion(problem)
res=steepestDescent(problem)

res$optValue
resY=(res$optY[[1]])
cutBy=svd(resY)$u[,1]
cutBy



##A is adjacency matrix and L is Laplacian matrix

par(mar=rep(0,4))
plot(c(1,5),c(1,(n/4+1)),type="n",bty="n",xaxt="n",yaxt="n")

for(i in 1:n){
  for(j in 1:n){
    if(A[i,j]==1){
      segments(px[i],py[i],px[j],py[j])
    }
  }
}

for(i in 1:n){
  if(cutBy[i]>0){
    points(px[i],py[i],pch=16,col="red",cex=3)
  }else{
    points(px[i],py[i],pch=16,col="green",cex=3)
  }
}


