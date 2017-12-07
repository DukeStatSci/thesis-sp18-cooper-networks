library(MASS)
library(MCMCpack)



#Specify example dataset
n=100
p=2
m=n-50
Z<-matrix(nrow=n,ncol=n)
Y<-matrix(nrow=n,ncol=n)
diag(Y)<--1
#for(i in seq(1,n,1)){
#  Y[i,-which(Y[i,]==-1)]<-sample(1:n-1,n-1,replace=F)
#}
#for(i in seq(1,n,1)){
#  for(j in seq(1,n,1)){
#    if(Y[i,j]>m){
#      Y[i,j]=0
#    }
#  }
#}
for(i in seq(1,n,1)){
  for(j in seq(1,n,1)){
    if(i!=j){
      Y[i,j]<-rpois(n=1,lambda=4)
    }
  }
}
X<-array(dim=c(n,n,p))
X[,,1]<-matrix(1,nrow=n,ncol=n)
diag(X[,,1])<--1
A<-matrix(1,nrow=n/2,ncol=n/2)
B<-matrix(0,nrow=n/2,ncol=n/2)
X[,,2]<-rbind(A,B)
diag(X[,,2])<--1
#Specify prior parameters
beta_0<-rep(0,p)
Sigma_0<-matrix(nrow=p,ncol=p)
Sigma_0[1,]<-c(2,1)
Sigma_0[2,]<-c(1,2)
v_0<-4
S_0<-matrix(nrow=2,ncol=2)
S_0[1,]<-c(2,1)
S_0[2,]<-c(1,2)
I2<-matrix(0,nrow=n-1,ncol=n-1)
diag(I2)<-1/6

#Set starting values
beta<-mvrnorm(n=1,mu=beta_0,Sigma=Sigma_0)
Sigma_ab<-riwish(v_0,solve(S_0))
a<-c()
b<-c()
for(i in seq(1,n,1)){
  ab<-mvrnorm(n=1,mu=c(0,0),Sigma=Sigma_ab)
  a<-c(a,ab[1])
  b<-c(b,ab[2])
}
Theta<-matrix(0,nrow=n,ncol=n)
diag(Theta)<--1
for(i in seq(1,n,1)){
  for(j in seq(1,n,1)){
    if(i!=j){
      Theta[i,j]<-beta%*%X[i,j,]+a[i]+b[j]
    }
  }
}
#Create matrices to store values
s=1000
betas<-matrix(nrow=s,ncol=p)
rs<-c()

#Start Gibbs Sampler
for(k in seq(1,s,1)){
  
  #Sample Sigma_ab from posterior
  Sigma_ab<-riwish(v_0+n,solve(S_0)+t(cbind(a,b))%*%cbind(a,b))
  sigma_a<-sqrt(Sigma_ab[1,1])
  sigma_b<-sqrt(Sigma_ab[2,2])
  
  #Sample a from posterior
  mu<-c()
  for(i in seq(1,n,1)){
    Theta2<-c()
    for(j in seq(1,n,1)){
      if(i!=j){
        Theta2<-c(Theta2,Theta[i,j]-t(beta)%*%X[i,j,]-b[j])
      }
    }
    mu<-c(mu,(sum(Theta2[-i])/(1/sigma_a^2+n-1)))
  }
  temp<-matrix(0,nrow=n,ncol=n)
  diag(temp)<-solve(1/sigma_a^2+n-1)
  a<-mvrnorm(n=1,mu=mu,Sigma=temp)
  
  #Sample b from posterior
  mu<-c()
  for(i in seq(1,n,1)){
    Theta2<-c()
    for(j in seq(1,n,1)){
      if(i!=j){
        Theta2<-c(Theta2,Theta[j,i]-t(beta)%*%X[j,i,]-a[i])
      }
    }
    mu<-c(mu,(sum(Theta2[-j])/(1/sigma_b^2+n-1)))
  }
  temp<-matrix(0,nrow=n,ncol=n)
  diag(temp)<-solve(1/sigma_b^2+n-1)
  a<-mvrnorm(n=1,mu=mu,Sigma=temp)
  
  
  #Sample beta from full conditional
  F1<-matrix(nrow=n^2-n,ncol=1)
  F1[,1]<-Theta[Theta!=-1]
  G<-matrix(nrow=n^2-n,ncol=2)
  H<-X[,,1]
  G[,1]<-H[H!=-1]
  H<-X[,,2]
  G[,2]<-H[H!=-1]
  beta_mu<-solve(solve(Sigma_0)+t(G)%*%G)%*%(solve(Sigma_0)%*%beta_0+t(G)%*%F1)
  beta_Sigma<-solve(solve(Sigma_0)+t(G)%*%G)
  
  beta<-mvrnorm(n=1,mu=beta_mu,Sigma=beta_Sigma)
  
  #Propose Thetas
  for(i in seq(1,n,1)){
    Xi<-X[i,,]
    Xi<-Xi[-i,]
    Yi<-Y[i,]
    Yi<-Yi[-i]
    Thetai<-Theta[i,]
    Thetai<-Thetai[-i]
    Tprop<-mvrnorm(n=1,mu=Xi%*%beta+a[i]+b[-i],Sigma=I2)
    logr<-sum(dpois(Yi,exp(Tprop),log=T))-
      sum(dpois(Yi,exp(Thetai),log=T))+
      sum(dnorm(Tprop,mean=Xi%*%beta+a[i]+b[-i],sd=sqrt(1/6),log=T))-
      sum(dnorm(Thetai,mean=Xi%*%beta+a[i]+b[-i],sd=sqrt(1/6),log=T))
    if(log(runif(1))<logr){
      Tprop<-append(Tprop,-1,after=i-1)
      Theta[i,]<-Tprop
      rs<-c(rs,1)
    }else{
      rs<-c(rs,0)
    }
  }
  
  #Store beta
  betas[k,]<-beta
}

print("Acceptance Rate")
print(mean(rs))
print("Posterior estimate for beta")
for(i in seq(1,p,1)){
  print(mean(betas[,i]))
}
print("Plots")
A<-seq(1,s,1)
plot(A,betas[,1],type="l", xlab="Iteration",ylab="Beta_1",main="Posterior Estimates from Gibbs Sampler for Beta_1")
plot(A,betas[,2],type="l", xlab="Iteration",ylab="Beta_2",main="Posterior Estimates from Gibbs Sampler for Beta_2")
acf(betas[,1],type="correlation",main="ACF for Beta_1")
acf(betas[,2],type="correlation",main="ACF for Beta_2")
