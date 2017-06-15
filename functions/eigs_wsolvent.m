function [KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, PHIP, CHI, IEIG)

RM = sqrt(r2(N));
Kmin = 1e-1/RM;
Kmax = 1e3/RM;

function EIGMIN = eig_wsolvent(K)
    [EIG,~]=gamma2_solvent(N,FA,K,CHI,PHIP);
    EIGMIN = EIG(IEIG);
end

KS = fminbnd(@eig_wsolvent,Kmin,Kmax);

[EIG,EIGV]=gamma2_solvent(N,FA,KS,CHI,PHIP);
EIGVS = EIGV(IEIG*2 - 1:IEIG*2);
EIGS = EIG(IEIG);

end
