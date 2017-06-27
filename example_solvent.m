clear;

%%%%%%%%% initializing %%%%%%%%%
addpath('functions')
addpath(genpath('chainstats'))
addpath('mkfigures')
addpath('misc')

%%%%%%%%% inputs %%%%%%%%%
% In the following inputs, Kuhn step b = 2*lp

% Number of Kuhn steps per polymer
N = 1;

% Flory-Huggins parameter between polymer and solvent
CHI = 0;

% Mole fraction of polymers
PHIP = 1e-10;

% Wavevector (q) in units of Kuhn steps, i.e. K = b*q
KV = logspace(-1, 5, 101);

%%%%%%%%% calculate structure factor %%%%%%%%%
S = scatter_homopoly(N, KV, CHI, PHIP);

%%%%%%%%% make a plot %%%%%%%%%
figure;hold
loglog(KV, S/PHIP)
xlabel('bq');ylabel('S(q)')
set(gca,'xscale','log');set(gca,'yscale','log');









%%% example of varying polymer concentration %%%
clear;
N = 5;
KV = logspace(-1, 2, 101);
CHI = 0;

PHIPV = [0.01,0.02,0.05,0.10];
figure;hold
for ii = 1:length(PHIPV)
    COL = (ii-1)/(length(PHIPV)-1);
    
    PHIP = PHIPV(ii)
    S = scatter_homopoly(N, KV, CHI, PHIP);
    
    plot(KV, S/PHIP, 'color', [COL 0 1-COL])
end
plotlog(1e1, 5e-1, -1, 1)  % add power law
xlabel('bq');ylabel('S(q)/\phi_P')
set(gca,'xscale','log');set(gca,'yscale','log');








clear;
% Take a look at Leo's data
file = load('data/Scattering Data.mat');
s = file.rawData;

inds = [19, 22, 1, 7, 10];
figure;hold
for ind = inds
    plot(s(:,ind), s(:,ind+1));
end
set(gca,'xscale','log');set(gca,'yscale','log')

% for ii = 1:8
%     ind1 = ii*3 - 2;
%     ind2 = ind1+1;
%     
%     figure;
%     loglog(s(:,ind1),s(:,ind2))
% end