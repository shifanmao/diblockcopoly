% test density-density correlations of diblock copolymers
addpath('../chainstats/')
addpath('../functions/')
addpath('../misc/')
addpath('../chainstats/eigcalc/')
clear

N=1e3;
FA=0.5;

% calculate density-density correlation at peak
[chis,ks,d2gam2]=spinodal(N,FA);

% calculate density-density at range of wavevectors
k=logspace(-1,2,100)*ks;
CHI=0.5*chis;
gam2=gamma2(N,FA,k,CHI);

% brazovskii estimation around peak
sq = 1./(-2*CHI+2*chis+0.5*d2gam2*power(k-ks,2));

% hohenberg-swift estimation around peak
eps4 = d2gam2./(8*k.^2);
sq2 = 1./(-2*CHI+2*chis+eps4.*power(k.^2-ks^2,2));

figure;
set(gca,'fontsize',18)
plot(k,1./gam2/N,k,sq/N,k,sq2/N,'linewidth',2)
set(gca,'xscale','log');set(gca,'yscale','log')
xlim([0.1,100]*ks);
% legend('S(q)/N','Brazovskii','Shankar')
legend('S(q)/N','Fredrickson-Helfand','Hohenberg-Swift')
xlabel('q');ylabel('S(q)/N')