N=0.56;
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
[chis,chi13,chi36,chi12,chi23,chi26]=plotphase(N,FAV);

figure;set(gca,'fontsize',18);hold
col1='k-';
col2='r-';
col3='b-';

% plot(FAV,chis*N,col1,'linewidth',2);
% plot(1-FAV,chis*N,col1,'linewidth',1.5)
% FA05 = find(abs(FAV-0.5)<1e-2);
% if ~isempty(FA05)
%     col = 'k';
%     plot(0.5,chis(FA05)*N,'o','color',col,...
%     'MarkerSize',8,'MarkerFaceColor',col)
% end

plot(FAV,(chi13+chis)*N,col2,'linewidth',1.5)
plot(1-FAV,(chi13+chis)*N,col2,'linewidth',1.5)

plot(FAV,(chi36+chis)*N,col3,'linewidth',1.5)
plot(1-FAV,(chi36+chis)*N,col3,'linewidth',1.5)

plot(FAV,(chi23(chi23>0)+chis(chi23>0))*N,col1,'linewidth',1.5)
plot(1-FAV,(chi23(chi23>0)+chis(chi23>0))*N,col1,'linewidth',1.5)


xlabel('f_A');ylabel('\chi N');box on
xlim([FAV(1),1-FAV(1)]);ylim([5,20]);