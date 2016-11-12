clear;
%close all
cd ../
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

% chain aspect ratio
alpha = 1;

% plot Ginzburg parameter
% Gaussian Chain limit
N = 1e6;
[chis,ks,d2gam2]=spinodal(N,0.5);
c = 0.5*d2gam2;

NV=logspace(-1,3,51);  % number of statistical steps of total chain
GiGCinv = ks^2*N/(4*pi)./power(alpha,3).*sqrt(1./NV/c);

% Wormlike Chain
NV=logspace(-1,3,51);  % number of statistical steps of total chain
[chis,ks,d2gam2]=spinodal(NV,0.5);
cV = 0.5*d2gam2;
GiWLCinv = ks.^2.*NV/(4*pi)./power(alpha,3).*sqrt(1./NV./cV);

%%%% Make a plot %%%%
figure;hold;
set(gca,'fontsize',20)
plot(NV,1./GiGCinv,'k-','linewidth',2)
plot(NV,1./GiWLCinv,'b-','linewidth',2)

% add power laws
x = logspace(2.4,3,10);y = power(x,1/2);
plot(x,0.45*y,'k--','linewidth',2)

x = logspace(-1,-0.5,10);y = power(x,2);
plot(x,0.1*y,'k--','linewidth',2)

set(gca,'xscale','log');set(gca,'yscale','log')
%set(gca,'ytick',[1e-2,1e-1,1e0,1e1,1e2,1e3])
xlabel('N');ylabel('Gi/\alpha^3');box on
cd mkfigures/