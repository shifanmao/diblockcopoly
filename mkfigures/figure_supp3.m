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

alphaV = [1,2,4,8,16];

% add asymptotes
for alpha = alphaV
    Nrange = logspace(-1, 0, 10);
    Yrange = 76.9*power(Nrange*alpha^(3/2), -4/3);
    plot(Nrange, Yrange, 'k--', 'linewidth', 3)
    
    Nrange = logspace(2, 3, 10);
    Yrange = 41.0*power(Nrange*alpha^6, -1/3);
    plot(Nrange, Yrange, 'k-.', 'linewidth', 3)    
end

% plot solutions
col = length(alphaV);
for jj = 1:length(alphaV)
    COL = (col-1)/(length(alphaV)-1);
    alpha = alphaV(jj);

    CHITN = fig1(fig1(:,1)==alpha, 3);
    plot(NV, CHITN-chis.*NV','-','color',[COL 0 1-COL],'linewidth',3);
    col = col-1;
end

% add asymptotes
for alpha = alphaV
    Nrange = logspace(-1, 0, 10);
    Yrange = 76.9*power(Nrange*alpha^(3/2), -4/3);
    plot(Nrange, Yrange, 'k--', 'linewidth', 3)
    
    Nrange = logspace(2, 3, 10);
    Yrange = 41.0*power(Nrange*alpha^6, -1/3);
    plot(Nrange, Yrange, 'k-.', 'linewidth', 3)    
end
legend('76.9(N\alpha^{3/2})^{-4/3}', '41(N\alpha^6)^{-1/3}')

xlabel('N');
ylabel('$\chi^{\mathrm{1L}}_{\mathrm{ODT}}N-\chi_{\mathrm{S}}^{\mathrm{MF}}N$',...
    'interpreter','latex')
set(gca,'xscale','log');set(gca,'yscale','log')
xlim([1e-1,1e3]);
set(gca,'linewidth',1.5)
box on

cd mkfigures/