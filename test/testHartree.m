% test validity of Hartree approximation
clear;close all
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

NV = logspace(0,4,10);
ht = zeros(length(NV),1);ii=1;
for N = NV
N
% CHI=1;
% tau=N*(2*CHI-chis);

FAV=0.5;
NQ=1;  % assume to Q dependence
[gam3,gam4]=calcgamma(N,FAV,NQ);

[chis,ks,d2gam2]=spinodal(N,FAV);
c = 0.5*d2gam2;
miu=N*gam3/power(c,3/2);
lam=N*gam4/power(c,2);

ht(ii) = power(lam*ks^2/4/pi,3/5)*power(ks,1/5);
ii = ii+1;
end

figure;hold
plot(NV,ht);

% add power law
x = logspace(3,5);
y = 32.024*power(x,-3/10);
plot(x,y,'k--');

cd test/