% calculate (test) vertex functions
% consider only magnitude dependence of gamma4(p,q)
%   and ignore angular dependence between p and q
clear all;close all
addpath('../chainstats/integrals/')
addpath('../chainstats/')
addpath('../functions/')
addpath('../misc/')

N = 1e4;
FA = 0.5;
NQ = 100;
CHI = 10/N;

% calculate spinodal
[~,ks]=spinodal(N,FA);
PMAG = linspace(0.5*ks,10*ks,NQ);

% calculate gamma4
gam4=zeros(1,NQ);
gam2=zeros(1,NQ);
for IQ = 1:NQ
    IQ
    K0 = [1,0,0];
    
    K1 = K0*ks;
    K2 = -K0*ks;
    K3 = K0*PMAG(IQ);
    K4 = -K0*PMAG(IQ);
    gam4(IQ) = gamma4(N,FA,K1,K2,K3,K4);
    gam2(IQ) = gamma2(N,FA,PMAG(IQ),CHI);
end