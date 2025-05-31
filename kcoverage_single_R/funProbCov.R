# funProbCov calculates and returns SINR-based coverage probability
# under arbtirary shadowing
# PCov=funProbCov(tValues,betaConst,x,numbMC,k)
# CovP is the 1-coverage probability
# tValues are the SINR threshold values. tValues can be a vector
# betaConst is the pathloss exponent.
# x is the input variable that incorporates the model parameters
# That is, x=W*a^(-2/betaConst)
# numbMc is number of sample (Sobol) points for  quasi-MC integration
# betaConst, x  and numbMC are scalars.
#
# Author: H.P. Keeler, Inria Paris/ENS, 2013
#
# References:
# [1] Keeler, B. Błaszczyszyn and M. Karray,
# 'SINR-based k-coverage probability in cellular networks with arbitrary
# shadowing', accepted at ISIT, 2013 
#
# R script originally based on Matlab script by H.P. Keeler, 2013
# hacked into R code by H.P. Keeler, 2013

funProbCov<-function(tValues,betaConst,x,numbMC,k){

  #General equation stems from Corollary 7 in [1]
  #source('funJn.R'); source('funIn.R');
  PCov<-rep(0,length(tValues));
  #maxium number connections possible; see Lemma in [1]
  nValues<-pmax(1,ceiling(1/tValues)); #note use of pmax
  for (j in 1:length(tValues)){
    nMax<-nValues[j];
    
    for (n in k:nMax){
      Tn<-tValues[j]/(1-tValues[j]*(n-1));     #eq. (17) in [1]   
      In<-funIn(betaConst,n,x); #eq. (12) in [1]
      Jn<-funJn(Tn,betaConst,n,numbMC); #eq. (15) in [1]
      Sn<-Tn^(-2*n/betaConst)*In*Jn; #eq. (18) in [1]
      PCov[j]<-(-1)^(n-k)*choose(n-1,k-1)*Sn+PCov[j]; #eq. (8) in [1]
    }   
  }
  funProbCov<-PCov;
}