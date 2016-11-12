clear;
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

%%%%%%%%%%%%%%%% SAME MOLECULAR WEIGHT EXAMPLES %%%%%%%%%%%%%%%%
N=10;  % number of statistical steps of total chain
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
alpha=4;  % dimensionless excluded volume parameter in the Gaussian limit
        % In the Gaussian chain limit, Nbar = C^2

% Figure 1: make a mean-field phase diagram
plotphase(N,FAV);

% Figure 2: make a phase diagram with density fluctuations
plotphaseRG(N,alpha,FAV);

% Figure 3-4: density-density correlations
densityRG(N,alpha,0.5);

%%%%%%%%%%%%%%%% MOLECULAR WEIGHT DEPENDENCE EXAMPLES %%%%%%%%%%%%%%%%
% Figure 5-6: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(-1,4,21)';  % number of statistical steps of total chain
chis=zeros(length(NV),1);
for ii = 1:length(NV)
    [chis(ii),ks,d2gam2]=spinodal(NV(ii),0.5);
end

figure;hold;set(gca,'fontsize',20);
plot(NV,chis.*NV,'linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');
xlabel('N');ylabel('\chiN');box on

figure;hold;set(gca,'fontsize',20);
plot(NV,1./ks,'linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','log');
xlabel('N');ylabel('1/q^*');box on

% Figure 7: renormalized spinodal at FA=0.5
chit=zeros(length(NV),1);
for ii = 1:length(NV)
    [chit(ii),phase]=spinodalRG(NV(ii),alpha,0.5);
end
figure;hold;set(gca,'fontsize',20)
col='b';
plot(NV,chis.*NV,'-','linewidth',2,'color',col)
plot(NV,chit.*NV,'s-','MarkerSize',8,'MarkerFaceColor',col,'MarkerEdgeColor',col);
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