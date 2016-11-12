clear;close all;
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

N=1e4;  % number of statistical steps of total chain
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
C=1e3;  % dimensionless excluded volume parameter in the Gaussian limit
        % In the Gaussian chain limit, Nbar = C^2

% Figure 3-4: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(-1,3,51);  % number of statistical steps of total chain
[chis,ks,d2gam2]=spinodal(NV,0.5);

figure;hold;set(gca,'fontsize',20);
plot(NV,chis.*NV,'linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');
xlabel('N');ylabel('\chiN');box on

% Figure 5: renormalized spinodal at FA=0.5
CV=logspace(1,4,21)';
[chit,phase]=spinodalRG(N,CV,0.5);
chit=reshape(chit,length(CV),1);

figure;hold;set(gca,'fontsize',20)
col='b';
plot(CV.^2,ones(length(CV),1)*spinodal(N,0.5)*N,'-','linewidth',2,'color',col)
plot(CV.^2,chit*N,'s','MarkerSize',8,'MarkerFaceColor',col,'MarkerEdgeColor',col);



figure;
plot(NV,chit*N-ones(length(CV),1)*spinodal(N,0.5)*N)