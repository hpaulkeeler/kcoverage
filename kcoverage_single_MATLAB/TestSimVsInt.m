% Calculates and plots SINR-based probability of k-coverage via simulation 
% of model and integration method outlined [1]
% 
%
% Author: H.P. Keeler, Inria Paris/ENS, 2014
%
% References:
% [1] H.P. Keeler, B. Blaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary 
% shadowing', ISIT, 2013
%
% NOTE: If you use this code in published work, please cite paper[1] by
% Keeler, Blaszczyszyn and Karray, as listed above.

clear all; close all; clc;

lambda=0.2887/2; %base station density

%K and betaConst values correspond to Walfisch-Ikegami model for a urban
%environment
betaConst=3.8; %path-loss exponent  
K=6910;

%log normal parameters
sigmDb=10;
sigma=sigmDb/10*log(10);
ESTwoBeta=exp(sigma^2*(2-betaConst)/betaConst^2);

%model constant - incorporates model parameters
a=lambda*pi*ESTwoBeta/K^2; %equation (6) in [1]

%noise parameters
N=10^(-96/10)/1000;
P=10^(62.2/10)/1000;
W=N/P;

%SINR threshold values 
tMinDb=-10;tMaxDb=0;
tValuesDb=(tMinDb:tMaxDb)'; %values in dB
tValues=10.^(tValuesDb/10);
tNumb=length(tValues);

%coverage number ie maximum number of connected base stations
k=3;
%analytic/integration section
numbMC=2^10; %number of integral sample points (based 2 for Sobol)
PCov=funProbCov(tValues,betaConst,W*a^(-betaConst/2),numbMC,k);

%Simulation section
simNumb=10^4; %number of simulations
diskRadius=20; %radius of simulation disk region (has to be larger when fading is incorporated)
PCovSim=funSimLogNormProbCov(tValues,betaConst,K,lambda,sigma,W,diskRadius,simNumb,k);

%plotting section

Pn=1-PCov;
PnSim=1-PCovSim;

%create suitable label
if W==0
    legendLabel='SIR';
else 
    legendLabel='SINR';
end

figure; 
plot(tValuesDb,PnSim,'o',tValuesDb,Pn);grid;
legend(['Simulation ', legendLabel],['Integration', legendLabel],'Location','NorthWest')
xlabel('T (dB)'); ylabel('1-P_c(T)')

