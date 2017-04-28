function [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP)

chis = 0/N/PHIP;
chie = 100/N/PHIP;
while (chie - chis) > 1e-4/N
    mid = (chie + chis) / 2;
    CHI = setCHIAB(mid);

    [~, MINEIG, ~] = eigs_wsolvent(N, FA, KV, PHIP, CHI, 1);
    if MINEIG*N < 1e-2
        chie = mid;
    else
        chis = mid;
    end
end
CHIABS = (chis + chie) / 2;

CHI = setCHIAB(CHIABS);
[KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, KV, PHIP, CHI, 1);

end

function CHI = setCHIAB(CHIAB)
    CHIBA = CHIAB;
    CHI = [0, CHIAB; CHIBA, 0];
end
% function CHIABS = spinodal_wsolvent(N,FA,KV,PHIP)
% 
% MINEIG = 1;
% CHIABS = 5/N;
% while MINEIG > 0
%     CHIABS = CHIABS + 0.01/N;
%     MINEIG = mineig(N, FA, KV, CHIABS, PHIP);
%     [N, CHIABS*N, MINEIG]
% end
% 
% end

function MINEIG = mineig(N, FA, KV, CHIAB, PHIP)
    CHIBA = CHIAB;
    CHI = [0, CHIAB; CHIBA, 0];

    [EIG,~,~]=gamma2_solvent(N,FA,KV,CHI,PHIP);
%     MINEIGS = min(EIG, [], 2);
    
    MINEIGS = EIG(:,1);
    [MINEIG, ~] = min(MINEIGS);
end

% function [CHIABS,KSS,EIGVS] = spinodal_wsolvent(N,FA,KV,PHIP)
% % find spinodal using binary search
% 
% chis = 0/N;
% chie = 500/N;
% while (chie - chis) > 1e-2/N
%     mid = (chie + chis) / 2;
%     CHIAB = mid;
%     CHIBA = CHIAB;
%     CHI = [0, CHIAB; CHIBA, 0];
%     
%     [EIG,~,~]=gamma2_solvent(N,FA,KV,CHI,PHIP);
%     MINEIGS = min(EIG, [], 2);
%     [MINEIG, KSIND] = min(MINEIGS);
%     
%     if MINEIG*N < 1e-5
%         chie = mid;
%     else
%         chis = mid;
%     end
% end
% CHIABS = (chis + chie) / 2;
% KSS = KV(KSIND);
% 
% [EIG,EIGV,~]=gamma2_solvent(N,FA,KSS,CHI,PHIP);
% [~, MININD] = min(EIG, [], 2);
% EIGVS = EIGV(MININD + (MININD - 1)*2 : MININD + (MININD - 1)*2 + 1);
% 
% end