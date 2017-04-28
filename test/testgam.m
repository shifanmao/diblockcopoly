clear;
%close all
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

NRR = 50;
FAV = 0.9;
NV = [1e-2,logspace(0,2,10),1e4];
NN = length(NV);

%NN = 1;
%NV = 1e-2;
chisv = zeros(NN,1);
ksv = zeros(NN,1);
d2gam2v = zeros(NN,1);
gam3v = zeros(NN,1);
gam4v = zeros(NN,1);
for I = 1:length(NV)
I	  
N = NV(I);

% Figure 8-9: vertex functions
NQ=1;  % number of wavevector sets in calculating GAM4
[chis,ks,d2gam2]=spinodal(N,0.5);
[gam3,gam4]=calcgamma(N,FAV,NQ);

if N<=1e-2
  ksv(I)=ks*NRR;
  d2gam2v(I)=d2gam2;
  chisv(I)=chis*NRR;
  %NRR*gam4/power(d2gam2/2,2)
  gam3v(I)=gam3*NRR
  gam4v(I)=gam4*NRR
else
  ksv(I)=ks*N;
  d2gam2v(I)=d2gam2;
  chisv(I)=chis*N;
  gam3v(I)=gam3*N
  gam4v(I)=gam4*N
end
    
end