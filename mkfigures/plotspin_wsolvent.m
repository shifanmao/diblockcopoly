function [FAV, PHIPV, CHISNPHIP, KSS, EIGVSS] = plotspin_wsolvent(N, LOAD)
    filename = sprintf('../data/SPIN_N%.2f.mat', N);

    if LOAD
        load(filename);
    else
%         FAV = linspace(.1, .9, 51);
%         PHIPV = linspace(.001, .999, 51);
        FAV = [0.3,0.5];
        PHIPV = linspace(0.1, 0.999, 20);
        CHISNPHIP = zeros(length(FAV), length(PHIPV));
        KSS = zeros(length(FAV), length(PHIPV));
        EIGVSS = zeros(length(FAV), length(PHIPV), 2);

        for jj = 1:length(PHIPV)
            PHIP = PHIPV(jj)
            for ii = 1:length(FAV)
                FA = FAV(ii);

                [CHIABS, KS, ~, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
                CHISNPHIP(ii,jj)=CHIABS*N*PHIP;
                KSS(ii,jj)=KS;
                EIGVSS(ii,jj,1:2)=EIGVS;
            end
        end
        save(filename, 'FAV', 'PHIPV', 'CHISNPHIP', 'KSS', 'EIGVSS'); 
    end
    
    %         figure;plotcontour(FAV, 1-PHIPV, CHISNPHIP)
%     figure;plotsurf(FAV, 1-PHIPV, CHISNPHIP)
%     figure;plotsurf(FAV, 1-PHIPV, KSS)
end