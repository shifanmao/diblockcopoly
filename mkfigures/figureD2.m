clear;
cd ..

alphaV = fliplr([1,2,3,4,5]);
NV=logspace(0,4,51)';

f = figure('Position', [0, 0, 1200, 800]);
hold;
set(gca,'fontsize',50)

% calculate mean-field results
chis = zeros(length(NV),1);
for ii = 1:length(NV)
    [chis(ii),ks,d2gam2]=spinodal(NV(ii),0.5);
end
plot(NV,chis.*NV,'k--','linewidth',3)

% % plot mean-field results
% NV1 = logspace(0,0.5,20);plot(NV1,ones(1,length(NV1))*6.135,'k.','linewidth',1.5)
% NV2 = logspace(2,3,20);plot(NV2,ones(1,length(NV2))*10.495,'k.','linewidth',1.5)
plot(1,6.135,'ko','markerfacecolor','k','markersize',8)
plot(500,10.495,'ko','markerfacecolor','k','markersize',8)

col = 1;
for jj = 1:length(alphaV)
    COL = (col-1)/(length(alphaV)-1);
    chit = zeros(1,length(NV));

    % plot one-loop results
    alpha = alphaV(jj);

    filename = strcat('data/chiODT',sprintf('_LP3overV%.2f',alpha^3));
    c1 = load(filename);
    plot(c1(:,1),c1(:,2),'-','color',[COL 0 1-COL],'linewidth',3);

    % plot Fredrickson-Helfand Theory
    chiNFH = 10.49+41.0*power(NV*alpha^6,-1/3);
    plot(NV,chiNFH,'--','color',[COL 0 1-COL],'linewidth',1.5)

    col = col+1;
end

xlabel('N');
ylabel('$\chi_{\mathrm{ODT}}^{\mathrm{1L}} N$','interpreter','latex')
ylabel('$(\chi N)_{\mathrm{ODT}}$','interpreter','latex')
set(gca,'xscale','log')
% set(gca,'ytick',[4,10,20,30,40,50])
set(gca,'ytick',[0,10,20,30,40,50])
set(gca,'xtick',[1,5,10,50,100,500])
xlim([1,5e2]);
ylim([0,30]);
set(gca,'linewidth',1.5)
box on

cd mkfigures/
