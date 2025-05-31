# funJn calculates and returns the In integral (eq (15) in [1])
# Jn=funJn(Tn,betaConst,n,numbMC)
# Jn is the J_{n,\beta} integral (scalar) value
# T_n is the n-based SINR threshold value (eq (17) in [1])
# betaConstant is path-loss exponent
# n is integer parameter
# betaConst, n, and numbMC are scalars. T_n can be a vector
# numbMc is number of sample (Sobol) points for  quasi-MC integration
#
# Author: H.P. Keeler, Inria Paris/ENS, 2013
#
# References
# [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
# 'SINR-based k-coverage probability in cellular networks with arbitrary
# shadowing', accepted at ISIT, 2013 
#
# R script originally based on Matlab script by H.P. Keeler, 2013
# hacked into R code by H.P. Keeler, 2013


funJn<-function(Tn,betaConst,n,numbMC){
#Calculates Jn with various integration methods
#n =2 and 3 uses quadrature methods respectively
#n>3 uses quasi Monte Carlo based on Sobol points
#function is called by funProbCov; see Corollary 7 in [1]

#Jn<-rep(0,length(Tn));
#for (k in 1:length(Tn)){    
  #Use  quadrature methods
  if (n==3){            
    fv3 <- function(arg, weight) {
      v1 <- arg[1];v2 <- arg[2];      
      ff <- (1/((v1*v2)+Tn))*(1/((v1*(1-v2))+Tn))*(v1*v2*(1-v2)*(1-v1))^(2/betaConst)*v1^(2/betaConst+1);
      return(ff)
    } # end integrand    
    Jn<-cuhre(2,1,fv3,lower=rep(0,2),upper=rep(1,2),flags= list(verbose=0, final=1))$value; #uses CUBA package
    
  } else if (n==2){     
    fv2 <- function(v1) (v1*(1-v1))^(2/betaConst)/(v1+Tn);        
    Jn <- integrate(fv2, lower=0, upper=1)$val #perform single quadrature
    
  }else if (n==1){
    Jn<-rep(1,length(Tn)); #return ones since J_1=1;
    
  }else {
    #Use QMC method        
    numbMCn<-(n-1)*numbMC; #scale number of points by dimension            
    # Use sobol points 
     qRandAll<-sobol(n = numbMCn, dim = n-1, scrambling = 3, normal = FALSE);
    # Use halton points
    #qRandAll<-halton(n = numbMCn, dim = n-1);
    
    eta_i<-vector('list',n);
    
    for (i in 1:n){
      if (i==1){
        eta_i[[i]]<-apply(qRandAll,1,prod);
      } else if (i==n){
        eta_i[[i]]<-(1-qRandAll[,n-1]);
      }else if (i==n-1){
        eta_i[[i]]<-(1-qRandAll[,(i-1)])*qRandAll[,i:(n-1)];
      }else{
        eta_i[[i]]<-(1-qRandAll[,(i-1)])*apply(qRandAll[,i:(n-1)],1,prod);
      } 
    }    
    
    # create/sample nominator and denominator of integral kernel
    numProdv_i<-rep(1,numbMCn);
    denomProdv_i<-numProdv_i;
    for (i in 1:(n-1)){
      viRand<-qRandAll[,i];      
      numProdv_i<-(viRand)^(i*(2/betaConst+1)-1)*(1-viRand)^(2/betaConst)*numProdv_i #numerator            
      denomProdv_i<-(eta_i[[i]]+Tn)*denomProdv_i; #denominator term
    }
    
    denomProdv_i=(eta_i[[n]]+Tn)*denomProdv_i;
    
    #factor out one term, arbitrarily choose the j-th term
    i=n-1;
    denomProdv_i=(eta_i[[i]]+Tn)/denomProdv_i;        
    
    kernelInt=numProdv_i*denomProdv_i; #integral kernerl
    Jn=mean(kernelInt); #perform (q)MC step
    
  }    
  
  funJn<-Jn;
}




