clear;close all
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

%%%%%%%%%%%%%%%% SAME MOLECULAR WEIGHT EXAMPLES %%%%%%%%%%%%%%%%
N=100;  % number of statistical steps of total chain
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
alpha=4;  % chain aspect ratio

% Figure 1: make a mean-field phase diagram
plotphase(N,FAV);

% Figure 2: make a phase diagram with density fluctuations
plotphaseRG(N,alpha,FAV);

% Figure 3-4: density-density correlations
densityRG(N,alpha,0.5);
figure(3);axis([0,1.5,0,.3])
figure(4);axis([9,12,0,3])

%%%%%%%%%%%%%%%% MOLECULAR WEIGHT DEPENDENCE EXAMPLES %%%%%%%%%%%%%%%%
% Figure 5-6: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(-1,4,21)';  % number of statistical steps of total chain
chis=zeros(length(NV),1);
ks=zeros(length(NV),1);
for ii = 1:length(NV)
    [chis(ii),ks(ii),d2gam2]=spinodal(NV(ii),0.5);
end

figure;hold;set(gca,'fontsize',20);
plot(NV,chis.*NV,'k-','linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');
xlim([NV(1),NV(end)])
xlabel('N');ylabel('\chiN');box on

figure;hold;set(gca,'fontsize',20);
plot(NV,1./ks,'k-','linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([NV(1),NV(end)])
xlabel('N');ylabel('1/q^*');box on

% Figure 7: renormalized ODTs at FA=0.5
figure;hold;set(gca,'fontsize',20)
plot(NV,chis.*NV,'k--','linewidth',2)

alpha=2;
chit=zeros(length(NV),1);
for ii = 1:length(NV)
    [chit(ii),phase]=spinodalRG(NV(ii),alpha,0.5);
end
plot(NV,chit.*NV,'r-','MarkerSize',8,'linewidth',2)

alpha=4;
chit=zeros(length(NV),1);
for ii = 1:length(NV)
    [chit(ii),phase]=spinodalRG(NV(ii),alpha,0.5);
end
plot(NV,chit.*NV,'b-','MarkerSize',8,'linewidth',2)
axis([1,500,4,30]);box on
set(gca,'xscale','log')

% Figure 8-9: vertex functions
NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
figure;hold;set(gca,'fontsize',20)
plot(FAV,-gam3*N,'k-','linewidth',2);xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)');box on

figure;hold;set(gca,'fontsize',20)
plot(FAV,gam4*N,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*)');box on

figure(1);saveas(gca,'example_figures/MFphase.png','png')
figure(2);saveas(gca,'example_figures/FLCphase.png','png')
figure(3);saveas(gca,'example_figures/psi2.png','png')
figure(4);saveas(gca,'example_figures/sinv.png','png')
figure(5);saveas(gca,'example_figures/spinodal.png','png')
figure(6);saveas(gca,'example_figures/domain.png','png')
figure(7);saveas(gca,'example_figures/FLCspinodal.png','png')
figure(8);saveas(gca,'example_figures/gam3.png','png')
figure(9);saveas(gca,'example_figures/gam4.png','png')