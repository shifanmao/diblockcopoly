% Plot density correlations at FA=0.5 (transition to LAM phase)
% Correspond to Figure 2 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"

clear
cd ..
addpath(genpath('.'))
NV = [100,10];
FA = 0.5;

%%%%%%%%%%%%%%%% alpha = 2 %%%%%%%%%%%%%%%%
alpha = 2;
p=[];

figure('position', [0,0,400,500]);
set(gca,'fontsize',15)
for ii = 1:length(NV)
    N = NV(ii);
    [k,Smf,Sfh,chit,chis,CHIV,sinvmf,sinvfh]=densityRG(N,alpha,FA,1,1);
    xlim([0,4])
end

% %%%%%%%%%%%%%%%% alpha = 4 %%%%%%%%%%%%%%%%
% alpha = 2;
% p=[];
% 
% for ii = 1:length(NV)
%     N = NV(ii);
% 
%     col = colv(ii);
%     [k,Smf,Sfh,chit,chis,CHIV,sinvmf,sinvfh]=densityRG(N,alpha,FA,0,1);
%     
%     plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col);
%     plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col);
% 
%     xlim([1,17]);ylim([0,20]);box on
%     xlabel('\chi N');
%     ylabel('$N\left<\tilde{\psi}(\vec{q}^*)\tilde{\psi}(-\vec{q}^*)\right>^{-1}$','Interpreter','latex')
% 
%     plot(chis*N,1./sinvmf,'o','color',col,...
%         'MarkerSize',10,'MarkerFaceColor',col);
%     plot(chit*N,1./sinvfh,'s','color',col,...
%     'MarkerSize',10,'MarkerFaceColor',col);
% %     p=[p,plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col)];
% %     p=[p,plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col)];
%     axis([8,15,0,5])
% end
% box on
% xlabel('\chi N');
% ylabel('$Nv \left< \tilde{\psi}(\vec{q}^*)\tilde{\psi}(-\vec{q}^*) \right>^{-1}$','Interpreter','latex')
% title('\alpha=2')
% 
cd mkfigures/