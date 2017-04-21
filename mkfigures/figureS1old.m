clear;
addpath('../functions/')
addpath(genpath('../chainstats'))

N = 1e1;
FA = 0.3;
KV = logspace(-2, 4, 100);

PHIP = 1-1e-4;
% % [CHIABS,KSS,EIGVS] = spinodal_wsolvent(N,FA,KV,PHIP);
% figure;hold
% CHIABS = spinodal_wsolvent(N,FA,KV,PHIP);

figure;hold
cnt = 1;
for CHIAB = [0, 50]/N
COL = (cnt -1) / (2-1)
CHIBA = CHIAB;
CHI = [0, CHIAB; CHIBA, 0];

[EIG,EIGV,KS]=gamma2_solvent(N,FA,KV,CHI,PHIP);

plot(KV, 1./EIG(:,1), '-','color', [COL 0 1-COL])
plot(KV, 1./EIG(:,2), '--','color', [COL 0 1-COL])
cnt = cnt + 1;
end

set(gca,'xscale','log')
set(gca,'yscale','log')

% CHI=0;
% val=gamma2(N,FA,kv,CHI);
% figure;
% plot(kv*N, 1./val)
% set(gca,'xscale','log')
% set(gca,'yscale','log')

% 
% s2invAA = [];
% s2AA = [];
% for k = kv
%     s2inv = s2inverse(N,FA,k);
%     s2 = s2gc(N,FA,k);
%     s2invAA = [s2invAA, s2inv(1, 1)];
%     s2AA = [s2AA, s2(1, 1)];
% end
% 
% figure;loglog(kv*N, s2invAA)
% figure;loglog(kv*N, s2AA)