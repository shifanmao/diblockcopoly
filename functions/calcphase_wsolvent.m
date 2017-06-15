function [FAVV, PHIPV, EIGV, EIG, KSV] = calcphase_wsolvent(N, CHI, IEIG, LOAD)

% filename = sprintf('../data/N%.2fCHIABN%.2fIEIG%d.mat', N, CHI(1,2)*N, IEIG);
filename = sprintf('../data/N%.2fCHIABNS%.2fIEIG%dTriMesh.mat', N, CHI, IEIG);

if LOAD
    load(filename);
else
    PHIPV = [0, linspace(.001, .999, 49), 1];
    FAV0 = .0;
    FAVf = 1;
    lenFAVV = zeros(length(PHIPV), 1);
    
    % triangular mesh
    FAVV = zeros(length(PHIPV), length(PHIPV));
    for jj = 1:length(PHIPV)
        if jj == 1  % PHIPV = 0
            lenFAVV(jj) = 1;
            FAVV(1, jj) = 0.5;
            FAVV(2:end, jj) = nan;
        else
            lenFAVV(jj) = jj;
            FAVV(1:lenFAVV(jj), jj) = linspace(FAV0, FAVf, lenFAVV(jj));
            FAVV(lenFAVV(jj)+1: end, jj) = nan;
        end
    end
    
    EIGV = zeros(length(PHIPV), length(PHIPV), 2);
    EIG = zeros(length(PHIPV), length(PHIPV));
    KSV = zeros(length(PHIPV), length(PHIPV));

    for jj = 1:length(PHIPV)
        PHIP = PHIPV(jj)
        
        for ii = 1:lenFAVV(jj)
            FA = FAVV(ii, jj);

            % use spinodal CHIAB
            if CHI == 0
                CHIAB = 0;
            else
                [CHIABS, ~, ~, ~] = spinodal_wsolvent(N, FA, PHIP);
                CHIAB = CHIABS*CHI;
            end
            CHIMAT = [0,CHIAB;CHIAB,0];

            [KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, PHIP, CHIMAT, IEIG);
            EIGV(ii, jj, 1:2) = EIGVS;
            EIG(ii, jj) = EIGS;
            KSV(ii, jj) = KS;
        end

        for ii = (lenFAVV(jj)+1):length(PHIPV)
            EIGV(ii, jj, 1:2) = [nan, nan];
            EIG(ii, jj) = nan;
            KSV(ii, jj) = nan;
        end
    end
    save(filename, 'FAVV', 'PHIPV', 'EIGV', 'EIG', 'KSV'); 
end


% function [FAV, PHIPV, EIGV, EIG, KSV] = calcphase_wsolvent(N, CHI, IEIG, LOAD)
% 
% % filename = sprintf('../data/N%.2fCHIABN%.2fIEIG%d.mat', N, CHI(1,2)*N, IEIG);
% filename = sprintf('../data/N%.2fCHIABNS%.2fIEIG%d.mat', N, CHI, IEIG);
% 
% if LOAD
%     load(filename);
% else
%     FAV = linspace(.1, .9, 41);
%     PHIPV = linspace(.001, .999, 51);
%     EIGV = zeros(length(FAV), length(PHIPV), 2);
%     EIG = zeros(length(FAV), length(PHIPV));
%     KSV = zeros(length(FAV), length(PHIPV));
% 
%     for jj = 1:length(PHIPV)
%         PHIP = PHIPV(jj)
%         for ii = 1:length(FAV)
%             FA = FAV(ii);
% 
%             % use spinodal CHIAB
%             if CHI == 0
%                 CHIAB = 0;
%             else
%                 [CHIABS, ~, ~, ~] = spinodal_wsolvent(N, FA, PHIP);
%                 CHIAB = CHIABS*CHI;
%             end
%             CHIMAT = [0,CHIAB;CHIAB,0];
% 
%             [KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, PHIP, CHIMAT, IEIG);
%             EIGV(ii, jj, 1:2) = EIGVS;
%             EIG(ii, jj) = EIGS;
%             KSV(ii, jj) = KS;
%         end
%     end
%     save(filename, 'FAV', 'PHIPV', 'EIGV', 'EIG', 'KSV'); 
% end