% Plot density correlations at FA=0.5 (transition to LAM phase)
% Correspond to Figure 2 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"

clear
cd ..
addpath(genpath('.'))
NV = [100,10];
colv = ['k','b'];

%%%%%%%%%%%%%%%% alpha = 2 %%%%%%%%%%%%%%%%
alpha = 4;
% Figure 1A
FA = 0.5;

colv = ['k','b'];
p=[];
for ii = 1:length(NV)
    N = NV(ii);

    col = colv(ii);
    [k,Smf,Sfh,chit,chis,CHIV,sinvmf,sinvfh]=densityRG(N,alpha,FA);
    p=[p,plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col)];
    p=[p,plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col)];
    axis([8,15,0,5])
end
box on
xlabel('\chi N');
ylabel('$Nv \left< \tilde{\psi}(\vec{q}^*)\tilde{\psi}(-\vec{q}^*) \right>^{-1}$','Interpreter','latex')

cd mkfigures/
