library(MASS)
library(MCMCpack)
library(truncnorm)
library(ggplot2)
Beta_0_true=-3.26
Beta_row_true=Beta_col=Beta_dyad1=Beta_dyad2=1
Sigma_ab_true=matrix(c(1,.5,.5,1),nrow=2,ncol=2)
Sigma_e_true=matrix(c(1,.9,.9,1),nrow=2,ncol=2)
s<-10000
num_sim<-8
Store<-data.frame(Est=1,Beta=-1,Type=factor("Poisson",level=c("Poisson","Binary")),Sim=4)
for(q in seq(1,num_sim,1)){
  #Simulate X
  n=100
  p=4
  X<-array(dim=c(n,n,p))
  X_row<-rnorm(n=n,mean=0,sd=1)
  X_col<-rnorm(n=n,mean=0,sd=1)
  for(i in seq(1,n,1)){
    X[i,-i,3]<-rep(X_row[i],n-1)
    X[-i,i,4]<-rep(X_col[i],n-1)
  }
  X[,,1]<-rnorm(n=n^2,mean=0,sd=1)
  Z<-matrix(nrow=n,ncol=n)
  Z[,]<-rbinom(n=n^2,size=1,prob=.5)
  X[,,2]<-Z[,]/.42
  diag(X[,,1])<--1
  diag(X[,,2])<--1
  diag(X[,,3])<--1
  diag(X[,,4])<--1
  
  Y<-matrix(nrow=n,ncol=n)
  for(i in seq(1,n,1)){
    for(j in seq(1,n,1)){
      if(i!=j){
        ab<-mvrnorm(mu=c(0,0),Sigma=Sigma_ab_true)
        ep<-mvrnorm(mu=c(0,0),Sigma=Sigma_e_true)
        Y[i,j]<-rpois(1,exp(sum(X[i,j,])+Beta_0_true+ab[1]+ab[2]))
      }
    }
  }
  S<-matrix(0,nrow=n,ncol=n)
  for(i in seq(1,n,1)){
    for(j in seq(1,n,1)){
      if(i!=j){
        if(Y[i,j]>0){
          S[i,j]<-1
        }
      }
    }
  }
  diag(S)<--1
  #Specify prior parameters
  beta_0<-rep(0,p)
  Sigma_0<-matrix(nrow=p,ncol=p)
  Sigma_0[1,]<-c(2,1,1,1)
  Sigma_0[2,]<-c(1,2,1,1)
  Sigma_0[3,]<-c(1,1,2,1)
  Sigma_0[4,]<-c(1,1,1,2)
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
  betas_poisson<-matrix(nrow=s,ncol=p)
  rs_poisson<-c()
  
  #Start Gibbs Sampler for Count Data
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
    G<-matrix(nrow=n^2-n,ncol=4)
    H<-X[,,1]
    G[,1]<-H[H!=-1]
    H<-X[,,2]
    G[,2]<-H[H!=-1]
    H<-X[,,3]
    G[,3]<-H[H!=-1]
    H<-X[,,4]
    G[,4]<-H[H!=-1]
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
        rs_poisson<-c(rs_poisson,1)
      }else{
        rs_poisson<-c(rs_poisson,0)
      }
    }
    
    #Store beta
    betas_poisson[k,]<-beta
  }
  
  #Specify prior parameters
  beta_0<-rep(0,p)
  Sigma_0<-matrix(nrow=p,ncol=p)
  Sigma_0[1,]<-c(2,1,1,1)
  Sigma_0[2,]<-c(1,2,1,1)
  Sigma_0[3,]<-c(1,1,2,1)
  Sigma_0[4,]<-c(1,1,1,2)
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
  betas_binary<-matrix(nrow=s,ncol=p)
  rs_binary<-c()
  Z<-S
  
  #Start Gibbs Sampler for Count Data
  for(k in seq(1,s,1)){
    
    #Sample Sigma_ab from posterior
    Sigma_ab<-riwish(v_0+n,solve(S_0)+t(cbind(a,b))%*%cbind(a,b))
    sigma_a<-sqrt(Sigma_ab[1,1])
    sigma_b<-sqrt(Sigma_ab[2,2])
    
    #Sample a from posterior
    mu<-c()
    for(i in seq(1,n,1)){
      Z2<-c()
      for(j in seq(1,n,1)){
        if(i!=j){
          Z2<-c(Z2,Z[i,j]-t(beta)%*%X[i,j,]-b[j])
        }
      }
      mu<-c(mu,(sum(Z2[-i])/(1/sigma_a^2+n-1)))
    }
    temp<-matrix(0,nrow=n,ncol=n)
    diag(temp)<-solve(1/sigma_a^2+n-1)
    a<-mvrnorm(n=1,mu=mu,Sigma=temp)
    
    #Sample b from posterior
    mu<-c()
    for(i in seq(1,n,1)){
      Z2<-c()
      for(j in seq(1,n,1)){
        if(i!=j){
          Z2<-c(Z2,Z[j,i]-t(beta)%*%X[j,i,]-a[i])
        }
      }
      mu<-c(mu,(sum(Z2[-j])/(1/sigma_b^2+n-1)))
    }
    temp<-matrix(0,nrow=n,ncol=n)
    diag(temp)<-solve(1/sigma_b^2+n-1)
    a<-mvrnorm(n=1,mu=mu,Sigma=temp)
    
    
    #Sample beta from full conditional
    F1<-matrix(nrow=n^2-n,ncol=1)
    F1[,1]<-Z[Z!=-1]
    G<-matrix(nrow=n^2-n,ncol=4)
    H<-X[,,1]
    G[,1]<-H[H!=-1]
    H<-X[,,2]
    G[,2]<-H[H!=-1]
    H<-X[,,3]
    G[,3]<-H[H!=-1]
    H<-X[,,4]
    G[,4]<-H[H!=-1]
    beta_mu<-solve(solve(Sigma_0)+t(G)%*%G)%*%(solve(Sigma_0)%*%beta_0+t(G)%*%F1)
    beta_Sigma<-solve(solve(Sigma_0)+t(G)%*%G)
    
    beta<-mvrnorm(n=1,mu=beta_mu,Sigma=beta_Sigma)
    
    for(i in seq(1,n,1)){
      for(j in seq(1,n,1)){
        if(i!=j){
          if(S[i,j]>0){
            Z[i,j]<-rtruncnorm(n=1,a=0,b=Inf,mean=X[i,j,]%*%beta+a[i]+b[j],sd=1)
          }else{
            Z[i,j]<-rtruncnorm(n=1,a=-Inf,b=0,mean=X[i,j,]%*%beta+a[i]+b[j],sd=1)
          }
        }
      }
    }
    #Store beta
    betas_binary[k,]<-beta
  }
  
  for(w in seq(1,p,1)){
    Store<-rbind(Store,c(mean(betas_poisson[,w]),w,"Poisson",q))
    Store<-rbind(Store,c(mean(betas_binary[,w]),w,"Binary",q))
  }
}

ggplot(data=Store[Store$Beta==1,], aes(x=Sim, y=as.numeric(Est), group=Type, color=Type)) +
  geom_point() + geom_hline(yintercept = 1) + labs(colour = "Likelihood", x="Simulation", y="Beta_1 Estimate", title="Simulated Posterior Estimates for Beta_1")

ggplot(data=Store[Store$Beta==2,], aes(x=Sim, y=as.numeric(Est), group=Type, color=Type)) +
  geom_line() +
  geom_point() + geom_hline(yintercept = 1) + labs(colour = "Likelihood", x="Simulation", y="Beta_2 Estimate", title="Simulated Posterior Estimates for Beta_2")

ggplot(data=Store[Store$Beta==3,], aes(x=Sim, y=as.numeric(Est), group=Type, color=Type)) +
  geom_line() +
  geom_point() + geom_hline(yintercept = 1) + labs(colour = "Likelihood", x="Simulation", y="Beta_3 Estimate", title="Simulated Posterior Estimates for Beta_3")

ggplot(data=Store[Store$Beta==4,], aes(x=Sim, y=as.numeric(Est), group=Type, color=Type)) +
  geom_line() +
  geom_point() + geom_hline(yintercept = 1) + labs(colour = "Likelihood", x="Simulation", y="Beta_4 Estimate", title="Simulated Posterior Estimates for Beta_4")


