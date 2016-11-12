% Plot density correlations at FA=0.5 (transition to LAM phase)
% Correspond to Figure 2 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"

cd ..
addpath('functions/')
fignum1 = 101;
fignum2 = 102;
NV = [100,10];
colv = ['k','b'];

%%%%%%%%%%%%%%%% alpha = 2 %%%%%%%%%%%%%%%%
fignum = fignum1;
LP3overV = 2^3;

% Figure 1A
FA = 0.5;
figure(fignum);
hold;set(gca,'fontsize',20)

colv = ['k','b'];
p=[];
for ii = 1:length(NV)
    N = NV(ii);
    C = power(sqrt(r2(N)),3)./N.*LP3overV;
    col = colv(ii);
%     [chis,chit,CHIV,Smf,Sfh,sinvmf,sinvfh]=densityRG(N,C,FA);
    [k,Smf,Sfh,chit,chis,CHIV,sinvmf,sinvfh]=densityRG(N,C,FA);
    figure(fignum);
    p=[p,plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col)];
    p=[p,plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col)];
    plot(chis*N,1./sinvmf,'o','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
    plot(chit*N,1./sinvfh,'s','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
end

figure(fignum);
%legend(p,{'N=10^3 Renormalized','N=10^3 Mean-field','N=1 Renormalized','N=1 Mean-field'})
title(strcat('\alpha',sprintf('=%d',power(LP3overV,1/3))))
axis([8,15,0,5])
box on
xlabel('\chi N');ylabel('$Nv \left< \tilde{\psi}(\vec{q}^*)\tilde{\psi}(-\vec{q}^*) \right>^{-1}$','Interpreter','latex')


%%%%%%%%%%%%%%%% alpha = 4 %%%%%%%%%%%%%%%%
fignum = fignum2;
LP3overV = 4^3;

% Figure 1A
FA = 0.5;
figure(fignum);
hold;set(gca,'fontsize',20)

p=[];
for ii = 1:length(NV)
    N = NV(ii);
    C = power(sqrt(r2(N)),3)./N.*LP3overV;
    col = colv(ii);
%     [chis,chit,CHIV,Smf,Sfh,sinvmf,sinvfh]=densityRG(N,C,FA);
    [k,Smf,Sfh,chit,chis,CHIV,sinvmf,sinvfh]=densityRG(N,C,FA);
    figure(fignum);
    p=[p,plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col)];
    p=[p,plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col)];
    plot(chis*N,1./sinvmf,'o','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
    plot(chit*N,1./sinvfh,'s','color',col,...
        'MarkerSize',8,'MarkerFaceColor',col);
end

figure(fignum);
%legend(p,{'N=10^3 Renormalized','N=10^3 Mean-field','N=1 Renormalized','N=1 Mean-field'})
title(strcat('\alpha',sprintf(' = %d',power(LP3overV,1/3))))
axis([8,15,0,5])
box on
xlabel('\chi N');ylabel('$Nv \left< \tilde{\psi}(\vec{q}^*)\tilde{\psi}(-\vec{q}^*) \right>^{-1}$','Interpreter','latex')





cd mkfigures/