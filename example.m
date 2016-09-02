clear;close all;
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

N=1e4;  % number of statistical steps of total chain
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
C=1e3;  % dimensionless excluded volume parameter in the Gaussian limit
        % In the Gaussian chain limit, Nbar = C^2

% Figure 1: make a mean-field phase diagram
plotphase(N,FAV);

% Figure 2: make a phase diagram with density fluctuations
plotphaseRG(N,C,FAV);

% Figure 3-4: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(-1,3,20);  % number of statistical steps of total chain
[chis,ks,d2gam2]=spinodal(NV,0.5);
figure;hold;set(gca,'fontsize',20);
plot(NV,chis.*NV,'linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');
xlabel('N');ylabel('\chiN')

figure;hold;set(gca,'fontsize',20);
plot(NV,1./ks,'linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','log');
xlabel('N');ylabel('1/q^*')

% Figure 5: renormalized spinodal at FA=0.5
CV=logspace(1,4,21)';
[chit,phase]=spinodalRG(N,CV,0.5);
chit=reshape(chit,length(CV),1);

figure;hold;set(gca,'fontsize',20)
col='b';
plot(CV.^2,ones(length(CV),1)*spinodal(N,0.5)*N,'-','linewidth',2,'color',col)
plot(CV.^2,chit*N,'s','MarkerSize',8,'MarkerFaceColor',col,'MarkerEdgeColor',col);

% %Empirical solutions
plot(CV.^2,ones(length(CV),1)*10.480+41.01*power(CV,-2/3),'-','linewidth',2,'color','k')
set(gca,'xscale','log');box on
xlabel('C^2');ylabel('\chiN');title(['N=',num2str(N)])
legend('MF theory','Renormalized ODT','Fit')

% Figure 6-7: density-density correlations
densityRG(N,C,0.5);

% Figure 8-9: vertex functions
NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
figure;hold;set(gca,'fontsize',20)
plot(FAV,-gam3*N,'k-','linewidth',2);xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)')

figure;hold;set(gca,'fontsize',20)
plot(FAV,gam4*N,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*)')

figure(1);saveas(gca,'example_figures/MFphase.png','png')
figure(2);saveas(gca,'example_figures/FLCphase.png','png')
figure(3);saveas(gca,'example_figures/spinodal.png','png')
figure(4);saveas(gca,'example_figures/domain.png','png')
figure(5);saveas(gca,'example_figures/FLCspinodal.png','png')
figure(6);saveas(gca,'example_figures/psi2.png','png')
figure(7);saveas(gca,'example_figures/sinv.png','png')
figure(8);saveas(gca,'example_figures/gam3.png','png')
figure(9);saveas(gca,'example_figures/gam4.png','png')