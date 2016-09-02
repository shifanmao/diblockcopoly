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

% Figure 3: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(-1,3,20);  % number of statistical steps of total chain
[chis,ks,d2gam2]=spinodal(NV,0.5);
figure;semilogx(NV,chis.*NV);xlabel('N');ylabel('\chiN')
figure;loglog(NV,1./ks);xlabel('N');ylabel('1/q^*')

% Figure 4: renormalized spinodal at FA=0.5
CV=logspace(1,4,21);
[chit,phase]=spinodalRG(N,CV,0.5);
chit=reshape(chit,length(CV),1);

figure;hold;set(gca,'fontsize',20)
col='b';
plot(CV.^2,ones(length(CV),1)*spinodal(N,0.5)*N,'--','linewidth',2,'color',col)
plot(CV.^2,chit*N,'s','MarkerSize',8,'MarkerFaceColor',col,'MarkerEdgeColor',col);

% %Empirical solutions
plot(CV.^2,ones(length(CV),1)*10.495,'--','linewidth',2,'color','k')
set(gca,'xscale','log');box on
xlabel('C^2');ylabel('\chiN');title(['N=',num2str(N)])
legend('MF theory','Renormalized ODT','Fit')

% Figure 5: density-density correlations
densityRG(N,C,0.5);

% Figure 6: vertex functions
NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4]=calcgamma(N,FAV,NQ);
figure;plot(FAV,-gam3*N,'k-','linewidth',2);xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)')
figure;plot(FAV,gam4*N,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*)')

% Save to images
saveas(gca,'../results/diblockcopoly-results/figures/phaseRGC1e2_flexible.eps','epsc')
