clear;
cd ..
addpath('functions/')
addpath(genpath('chainstats/'))
addpath('misc/')

alphaV = fliplr([1,2,3,4,5]);
NV=logspace(0,4,51)';

% f = figure('Position', [0, 0, 1200, 800]);
figure;
hold;
% set(gca,'fontsize',50)
set(gca,'fontsize',18)

% calculate mean-field results
chis = zeros(length(NV),1);
for ii = 1:length(NV)
    [chis(ii),ks,d2gam2]=spinodal(NV(ii),0.5);
end
plot(NV,chis.*NV,'k--','linewidth',3)

col = 1;
for jj = 1:length(alphaV)
    COL = (col-1)/(length(alphaV)-1);
    chit = zeros(1,length(NV));

    % plot one-loop results
    alpha = alphaV(jj);

    filename = strcat('data/chiODT',sprintf('_LP3overV%.2f',alpha^3));
    c1 = load(filename);
    plot(c1(:,1),c1(:,2),'-','color',[COL 0 1-COL],'linewidth',3);
    col = col+1;
end

% % plot mean-field results
% NV1 = logspace(0,0.5,20);plot(NV1,ones(1,length(NV1))*6.135,'k.','linewidth',1.5)
% NV2 = logspace(2,3,20);plot(NV2,ones(1,length(NV2))*10.495,'k.','linewidth',1.5)
plot(1,6.135,'ko','markerfacecolor','k','markersize',8)
plot(1000,10.495,'ko','markerfacecolor','k','markersize',8)

xlabel('N');
ylabel('$\chi^{\mathrm{1L}}_{\mathrm{ODT}}N$','interpreter','latex')
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xtick',[1,5,10,50,100,500,1000])
set(gca,'ytick', 5:5:30)
xlim([1,1e3]);
ylim([5,30]);
set(gca,'linewidth',1.5)
box on

% add asymptotes
col = 1;
for jj = 1:length(alphaV)
    COL = (col-1)/(length(alphaV)-1);
    
    alpha = alphaV(jj);
    
    Nrange = logspace(log10(50), 3, 10);
    Yrange = 10.495 + 41.0*power(Nrange*alpha^6, -1/3);
    plot(Nrange, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
    
%     Nrange = logspace(0, log10(3), 10);
    Nrange = logspace(0, log10(3), 10)/alpha^2;
    Yrange = 6.135 + 128.1*power(Nrange*alpha^2, -1);
    plot(Nrange, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
%     plot(Nrange, Yrange, 'k--')
    col = col+1;
end

col = 1;
for jj = 1:length(alphaV)
    COL = (col-1)/(length(alphaV)-1);    
    alpha = alphaV(jj);
    
%     Nrange = logspace(0, log10(5), 10)/alpha^2;
    Nrange = logspace(0, log10(5), 10);
    Yrange = 6.135 + 128.1*power(Nrange*alpha^2, -1);
    plot(Nrange, Yrange, '--', 'color', [COL 0 1-COL], 'linewidth', 3)
    col = col+1;
end


cd mkfigures/
