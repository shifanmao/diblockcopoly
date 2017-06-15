function [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP)

chis = 0/N/PHIP;
chie = 1e2/N/PHIP;

while (chie - chis) > 1e-2/N/PHIP
    mid = (chie + chis) / 2;
    CHI = setCHIAB(mid);

    [~, MINEIG, ~] = eigs_wsolvent(N, FA, PHIP, CHI, 1);
    if MINEIG*N < 1e-2
        chie = mid;
    else
        chis = mid;
    end
end
CHIABS = (chis + chie) / 2;
CHIABS = chis;

CHI = setCHIAB(CHIABS);
[KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, PHIP, CHI, 1);

end

function CHI = setCHIAB(CHIAB)
    CHIBA = CHIAB;
    CHI = [0, CHIAB; CHIBA, 0];
end