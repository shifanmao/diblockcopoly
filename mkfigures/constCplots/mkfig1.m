% Plot density correlations at FA=0.5 (transition to LAM phase)
% Correspond to Figure 2 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"

cd ..
addpath('functions/')
fignum1 = 101;
fignum2 = 102;
C1 = 1e3;
C2 = 1e4;

% Figure 1A
FA = 0.5;
C = sqrt(C1);
figure(fignum1);
hold;set(gca,'fontsize',20)

NV = [1e3,1];
colv = ['k','b'];
p=[];
for ii = 1:length(NV)
    N = NV(ii);
    col = colv(ii);
    [chis,chit,CHIV,Smf,Sfh,sinvmf,sinvfh]=densityRG(N,C,FA);
    figure(fignum1);
    p=[p,plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col)];
    p=[p,plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col)];
    plot(chis*N,1./sinvmf,'o','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
    plot(chit*N,1./sinvfh,'s','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
end

figure(fignum1);
legend(p,{'N=10^3 Renormalized','N=10^3 Mean-field','N=1 Renormalized','N=1 Mean-field'})
xlim([1,17]);ylim([0,25]);box on
xlabel('\chi N');ylabel('$N<\tilde{\psi}(q^*)\tilde{\psi}(-q^*)>^{-1}$','Interpreter','latex')
title(strcat('C^2=10^',sprintf('%d',log10(C1))))
savename = sprintf('mkfigures/figure1A.eps');
saveas(gcf,savename,'epsc')

% Figure 1B
FA = 0.5;
C = sqrt(C2);
figure(fignum2);
hold;set(gca,'fontsize',20)

NV = [1e3,1];
colv = ['k','b'];
p=[];
for ii = 1:length(NV)
    N = NV(ii);
    col = colv(ii);
    [chis,chit,CHIV,Smf,Sfh,sinvmf,sinvfh]=densityRG(N,C,FA);
    figure(fignum2);
    p=[p,plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col)];
    p=[p,plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col)];
    plot(chis*N,1./sinvmf,'o','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
    plot(chit*N,1./sinvfh,'s','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
end

figure(fignum2);
legend(p,{'N=10^3 Renormalized','N=10^3 Mean-field','N=1 Renormalized','N=1 Mean-field'})
xlim([1,17]);ylim([0,25]);box on
xlabel('\chi N');ylabel('$N<\tilde{\psi}(q^*)\tilde{\psi}(-q^*)>^{-1}$','Interpreter','latex')
title(strcat('C^2=10^',sprintf('%d',log10(C2))))
savename = sprintf('mkfigures/figure1B.eps');
saveas(gcf,savename,'epsc')

cd mkfigures