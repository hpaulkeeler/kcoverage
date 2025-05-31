# Calculates and plots SINR-based probability of k-coverage via simulation 
# of model and integration method outlined [1]
# 
#
# Author: H.P Keeler, INRIA Rocquencourt/ENS, 2013
#
# References
# [1] H.P Keeler, B. Błaszczyszyn and M. Karray,
# 'SINR-based coverage probability in cellular networks under multiple
#  connections', submitted to ISIT, 2013 
#
# R script originally based on Matlab script by H.P. Keeler, 2013
# hacked into R code by H.P. Keeler, 2013


rm(list=ls(all=TRUE));

library("R2Cuba"); #for integration  "cuhre" in funJn.R
library("randtoolbox") ; #for integration  "sobol" in funJn.R
codeDir<-"/home/ROCQ/trec/keeler/Documents/R/";
setwd(codeDir); #set source path
source(paste(codeDir,"funProbCov.R",sep=""));
source(paste(codeDir,"funIn.R",sep=""));
source(paste(codeDir,"funJn.R",sep=""));
source("funSimLogNormProbCov.R");


lambda<-0.2887/2; #base station density
lambda<-4.6188;
#K and betaConst values correspond to Walfisch-Ikegami model for a urban
#environment
betaConst<-3.8; #path-loss exponent  
K<-6910;

#log normal parameters
sigmDb<-10;
sigma<-sigmDb/10*log(10);
ESTwoBeta<-exp(sigma^2*(2-betaConst)/betaConst^2);

#model constant - incorporates model parameters
a<-lambda*pi*ESTwoBeta/K^2; #equation q (6) in [1]

#noise paramters
N<-10^(-96/10)/1000;
P<-10^(62.2/10)/1000;
W<-0;N/P;

#SINR threshold values 
tMinDb<--10;tMaxDb<-25;
tValuesDb<-(tMinDb:tMaxDb); #values in dB
tValues<-10.^(tValuesDb/10);
tNumb<-length(tValues);

#coverage number
k<-1;
#analytic/integration section
numbMC<-10^3;
PCov<-funProbCov(tValues,betaConst,W*a^(-betaConst/2),numbMC,k);

#Simulation section
simNumb<-10^4; #number of simulations
diskRadius<-20; #radius of simulation disk region (has to be larger when fading is incorporated)
PCovSim<-funSimLogNormProbCov(tValues,betaConst,K,lambda,sigma,W,diskRadius,simNumb,k);

#plotting section

Pn<-1-PCov;
PnSim<-1-PCovSim;

#pdf('SINR.pdf')
plot(tValuesDb,Pn,type='l',xlim=c(tValuesDb[1],tValuesDb[tNumb]),ylim=c(0,1),xlab='T (dB)',ylab='1-P_c(T)');
points(tValuesDb,PnSim);
#grid();
#dev.off();
