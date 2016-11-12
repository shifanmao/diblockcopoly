clear;
cd ..

PV = fliplr([1,2,3,4,5].^3);

f = figure('Position', [0, 0, 1200, 800]);
hold;
set(gca,'fontsize',50)

% calculate mean-field results
LP3overV = PV(1);
filename = strcat('data/chiODT',sprintf('_LP3overV%.2f',LP3overV));
c1 = load(filename);
NV = c1(:,1);
[chis,ks,d2gam2]=spinodal(NV,0.5);

col = 1;
for jj = 1:length(PV)
    COL = (col-1)/(length(PV)-1);
    chit = zeros(1,length(NV));

    % plot one-loop results
    LP3overV = PV(jj)
    alpha = power(LP3overV,1/3);
    
    filename = strcat('data/chiODT',sprintf('_LP3overV%.2f',LP3overV));
    c1 = load(filename);
    NV = c1(:,1);
    [chis,ks,d2gam2]=spinodal(NV,0.5);

    plot(c1(:,1),c1(:,2)-chis'.*c1(:,1),'-','color',[COL 0 1-COL],'linewidth',3);
    % add power laws
    x = logspace(1.5,3,10);y = power(x*power(alpha,6),-1/3);
    plot(x,41*y,'--','color',[COL 0 1-COL],'linewidth',3) % Fredrickson-Helfand Theory
    x = logspace(0,0.5,10);y = power(x*power(alpha,2),-1);
    plot(x,128*y,'--','color',[COL 0 1-COL],'linewidth',3)
    
%     plot(c1(:,1)*alpha^6,c1(:,2)-chis'.*c1(:,1),'-','color',[COL 0 1-COL],'linewidth',3);
    % add power laws
%     x1 = 1/alpha^6;
%     x = logspace(log10(x1),0.5,10);y = power(x*power(alpha,2),-1);
%     plot(x*alpha^6,128*y,'--','color',[COL 0 1-COL],'linewidth',3)

    col = col+1;
end

xlabel('N');
ylabel('$\chi_{\mathrm{ODT}}^{\mathrm{1L}}N-\chi_{\mathrm{S}}^{\mathrm{MF}}N$','interpreter','latex')
set(gca,'xscale','log')
set(gca,'yscale','log')

set(gca,'ytick',[1e-1,1e0,1e1,1e2])
% set(gca,'xtick',[1,5,10,50,100,500])
xlim([1,1e6])
ylim([1e-1,2e2])
set(gca,'linewidth',1.5)
box on

cd mkfigures/