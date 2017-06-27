clear;
cd ..
addpath('functions/')
addpath(genpath('chainstats/'))
addpath('misc/')

NV = logspace(-1,3,81);

figure;
hold;
set(gca,'fontsize',18)

% calculate mean-field results
chis = zeros(length(NV),1);
for ii = 1:length(NV)
    [chis(ii),ks,d2gam2]=spinodal(NV(ii),0.5);
end

% plot one-loop results
filename = 'fig1plot';
fig1 = load(filename);

% alphaV = fliplr([1,2,3,4,5,6]);
alphaV = fliplr([1,2,4,8,16]);
col = 1;
for jj = 1:length(alphaV)
    COL = (col-1)/(length(alphaV)-1);
    alpha = alphaV(jj);

    CHITN = fig1(fig1(:,1)==alpha, 3);
    plot(NV, CHITN,'-','color',[COL 0 1-COL],'linewidth',3);
%     plot(NV*alpha^2, CHITN,'-','color',[COL 0 1-COL],'linewidth',3);
    col = col+1;
end

% % plot mean-field results
plot(NV,chis.*NV','k--','linewidth',3)
plot(1e-1,6.135,'ko','markerfacecolor','k','markersize',8)
plot(1e3,10.495,'ko','markerfacecolor','k','markersize',8)

xlabel('N');
ylabel('$\chi^{\mathrm{1L}}_{\mathrm{ODT}}N$','interpreter','latex')
set(gca,'xscale','log');
set(gca,'yscale','linear');
% set(gca,'xtick',[1,5,10,50,100,500,1000])
set(gca,'ytick', 5:5:30)
xlim([1e-1,1e3]);
ylim([5,30]);
set(gca,'linewidth',1.5)
box on
% 
% % add asymptotes
% col = 1;
% for jj = 1:length(alphaV)
%     COL = (col-1)/(length(alphaV)-1);
%     
%     alpha = alphaV(jj);
%     
%     Nrange = logspace(log10(50), 3, 10);
%     Yrange = 10.495 + 41.0*power(Nrange*alpha^6, -1/3);
%     plot(Nrange, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
%     
% %     Nrange = logspace(0, log10(3), 10);
%     Nrange = logspace(0, log10(50), 10)/alpha^2;
%     Yrange = 6.135 + 111*power(Nrange*alpha^2, -1);
% %     plot(Nrange, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
%     plot(Nrange*alpha^2, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
%     col = col+1;
% end
% 
% col = 1;
% for jj = 1:length(alphaV)
%     COL = (col-1)/(length(alphaV)-1);    
%     alpha = alphaV(jj);
%     
% %     Nrange = logspace(0, log10(5), 10)/alpha^2;
%     Nrange = logspace(0, log10(5), 10);
%     Yrange = 6.135 + 128.1*power(Nrange*alpha^2, -1);
%     plot(Nrange, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
%     col = col+1;
% end


cd mkfigures/