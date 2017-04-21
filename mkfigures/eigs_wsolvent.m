function [KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, KV, PHIP, CHI, IEIG)

% [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);
% semilogx(KV, 1./EIG(:,IEIG))/N

% % IEIG = 1;
% MININD = find(1./EIG(:,IEIG)==max(1./EIG(:,IEIG)));
% KS = KV(MININD);
% EIGVS = EIGV(MININD, IEIG*2 - 1:IEIG*2);
% EIGS = EIG(MININD, IEIG);

Kmin = KV(1);
Kmax = KV(end);

function EIGMIN = eig_wsolvent(K)
    [EIG,~]=gamma2_solvent(N,FA,K,CHI,PHIP);
    EIGMIN = EIG(IEIG);
end

KS = fminbnd(@eig_wsolvent,Kmin,Kmax);

[EIG,EIGV]=gamma2_solvent(N,FA,KS,CHI,PHIP);
EIGVS = EIGV(IEIG*2 - 1:IEIG*2);
EIGS = EIG(IEIG);

end

