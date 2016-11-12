clear;close all;
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

N = 1e5;
FA = 0.5;

alphaV = logspace(-0.5,0.5,100)';
Nbar = N*power(alphaV,6);
[chit,phase]=spinodalRG(N,alphaV,0.5);
chit=reshape(chit,length(alphaV),1);

% analytical
FH = 10.495+41.0*power(Nbar,-1/3);

figure;
hold;
plot(Nbar,chit*N,'k.')
plot(Nbar,FH,'k-')
set(gca,'xscale','log')
% 
% % scattering functions
% alpha = 4;
% N = 10;
% densityRG(N,alpha,FA)