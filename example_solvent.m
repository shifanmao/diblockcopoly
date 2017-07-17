clear;
% 
% %%%%%%%%% initializing %%%%%%%%%
% addpath('functions')
% addpath(genpath('chainstats'))
% addpath('mkfigures')
% addpath('misc')
% 
% %%%%%%%%% inputs %%%%%%%%%
% % In the following inputs, Kuhn step b = 2*lp
% 
% figure;hold
% PHIPV = [1, 2, 3]*1e-1;
% 
% % Number of Kuhn steps per polymer
% N = 1;
% 
% % Flory-Huggins parameter between polymer and solvent
% CHI = .5/N/PHIP;
% 
% cnt = 1;
% for PHIP = PHIPV;
% col = (cnt-1)/(length(PHIPV)-1);
% % % Mole fraction of polymers
% % PHIP = 1e-10;
% 
% % Wavevector (q) in units of Kuhn steps, i.e. K = b*q
% KV = logspace(-1, 5, 101);
% 
% %%%%%%%%% calculate structure factor %%%%%%%%%
% S = scatter_homopoly(N, KV, CHI, PHIP);
% 
% %%%%%%%%% make a plot %%%%%%%%%
% 
% loglog(KV, S/PHIP,'color', [col 0 1-col]);
% 
% cnt = cnt+1;
% end
% 
% xlabel('bq');ylabel('S(q)')
% set(gca,'xscale','log');set(gca,'yscale','log');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%% example of varying polymer concentration %%%
% clear;
% N = 5;
% KV = logspace(-1, 2, 101);
% CHI = 0;
% 
% PHIPV = [0.01,0.02,0.05,0.10];
% figure;hold
% for ii = 1:length(PHIPV)
%     COL = (ii-1)/(length(PHIPV)-1);
%     
%     PHIP = PHIPV(ii)
%     S = scatter_homopoly(N, KV, CHI, PHIP);
%     
%     plot(KV, S/PHIP, 'color', [COL 0 1-COL])
% end
% plotlog(1e1, 5e-1, -1, 1)  % add power law
% xlabel('bq');ylabel('S(q)/\phi_P')
% set(gca,'xscale','log');set(gca,'yscale','log');


clear;
file = load('data/Scattering Data v2.mat');
s = file.rawData;
sl = file.rawDataLabels;

s2 = file.rawData2;
sl2 = file.rawData2Labels;


for ii = 1:length(sl)
    disp(strcat(num2str(ii), ' -> ', sl{ii}))
end

for ii = 1:length(sl2)
    disp(strcat(num2str(ii), ' -> ', sl2{ii}))
end






%%%%%%%%%% SAMPLE 1 %%%%%%%%%%%
inds = [19, 22, 1, 7, 10];
conc = [1, 2, 4, 7, 11];
s0 = 6e2;
L = 300;
lp = 300*5/10; % unit in Angstrom
N=L/(2*lp);

figure;hold;cnt = 1;
for ind = inds
    col = (cnt-1)/(length(inds)-1);
    plot(s(:,ind), s(:,ind+1)/conc(cnt), ...
        'color', [col 0 1-col],'linewidth',2);
    cnt = cnt+1;
end
set(gca,'xscale','log');set(gca,'yscale','log')

PHIPV = 8e-2*[conc];
PHIPV = 8e-2*[1,1];
KV = logspace(-1, 3, 101);
cnt = 1;
for PHIP = PHIPV
    col = (cnt-1)/(length(inds)-1);
    CHI=0/PHIP/N;

    S = scatter_homopoly(N, KV, CHI, PHIP);
    plot(KV/(2*lp),S/PHIP*s0/N,'--',...
        'color', [col 0 1-col],'linewidth',2);
    cnt = cnt+1;
end

plotlog(4e-2,2e2,-1,.5)

axis([1e-3,1,1e-1,1e4])
box on
set(gca,'fontsize',20)
xlabel('q(1/A)');ylabel('S(q)')
title('p-xylene')
legend('1mg/mL','2mg/mL','4mg/mL','7mg/mL','11mg/mL', 'location', 'southwest')













%%%%%%%%%% SAMPLE 2 %%%%%%%%%%%

inds = fliplr([1, 4, 7, 10, 13]);
conc = fliplr([11, 7, 4, 2, 1]);
s0 = 3e3;
L = 300;
% lp = 150; % unit in Angstrom
lp = 100; % unit in Angstrom
N=L/(2*lp);

figure;hold;cnt = 1;
for ind = inds
    col = (cnt-1)/(length(inds)-1);
    plot(s(:,ind), s(:,ind+1)/conc(cnt), ...
        'color', [col 0 1-col],'linewidth',2);
    cnt = cnt+1;
end
set(gca,'xscale','log');set(gca,'yscale','log')

PHIPV = 8e-2*[conc];
PHIPV = 8e-2*[1,1];
KV = logspace(-1, 3, 101);
cnt = 1;
for PHIP = PHIPV
    col = (cnt-1)/(length(inds)-1);
    CHI=0/PHIP/N;

    S = scatter_homopoly(N, KV, CHI, PHIP);
    plot(KV/(2*lp),S/PHIP*s0/N,'--',...
        'color', [col 0 1-col],'linewidth',2);
    cnt = cnt+1;
end

plotlog(4e-2,1e3,-1,.5)

axis([1e-3,1,1e-1,1e4])
box on
set(gca,'fontsize',20)
xlabel('q(1/A)');ylabel('S(q)')
title('Tetralin')
legend('1mg/mL','2mg/mL','4mg/mL','7mg/mL','11mg/mL', 'location', 'southwest')
