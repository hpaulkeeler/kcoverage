
# funSimLogNormProbCov returns SINR-based k-coverage probability 
# under log-normal shadowing based on repeated simulations of
# of model outlined in [1]
#
# simPCovk=funSimLogNormProbCov(tValues,betaConst,K,lambda,sigma,W,diskRadius,simNumb,k)
# simPCovk is the k-coverage probability
# tValues are the SINR threshold values. tValues can be a vector
# betaConst is the pathloss exponent.
# K = path-loss constant
# lambda = density of base station/nodes of cellular network
# sigma = standard deviation of log-normal shadowing (not in dB)
# W = noise constant
# diskRadius = radius of simulation region (Warning: if too small, edge
# effects will not be sufficiently included and disagreement with analytic
# results may occurr)
# simNumb = number of simulations
# k = coverage number (often set to one)
# All input values (except tValues) are scalars
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

funSimLogNormProbCov<-function(tValues,betaConst,K,lambda,sigma,W,diskRadius,simNumb,k){

#tValues=2;betaConst=4;K=3000;lambda=3;sigma=1/10;W=0;diskRadius=20; simNumb=10^3;k=1;

tNumb<-length(tValues);

### Simulation Section ###
#(uniformally) randomly places nodes on a disk of radius diskRadius
diskArea<-pi*diskRadius^2;
coveredNumbk<-rep(0,tNumb);
ESTwoBeta<-exp(sigma^2*(2-betaConst)/betaConst^2);
#rescale lambda - see foot note 5 in [1]
lambdaSim<-lambda*ESTwoBeta;
for (i in 1:simNumb){
  randNumb<-rpois(1,lambdaSim*diskArea); #poisson random variable
  #shadowing distribution can be constant if lambda is rescaled - see [1]
  shadowRand<-rep(1,randNumb); 
  #random distances from the typical node 
  rRand<-diskRadius*sqrt(runif(randNumb)); #uniform in cartesion, not polar coordinates
  
  signalRand<-shadowRand*(K*rRand)^(-betaConst);
  interferTotal<-sum(signalRand); #total inteference in network
  SINR<-signalRand/((interferTotal-signalRand)+W); #calculate SINR for each node in the network
  
  for (j in 1:tNumb){
    T<-tValues[j];
    #counts how many nodes are exactly k or more connected/covered
    if (sum(SINR>=T)>=k){
      coveredNumbk[j]<-coveredNumbk[j]+1;
    }
  }
  
}
simPCovk<-coveredNumbk/simNumb;


funSimLogNormProbCov<-simPCovk
}