% Figure 1: Spinodals of diblock copolymer solutions
clear;

CALCON = 0;
filename1 = '../data/figureS1A';
filename2 = '../data/figureS1B';

if CALCON

    FAV = linspace(.01, .99, 101);
    
%     % Figure 1A
%     N = 1e3;
%     PHIP = 0.5;
% 
%     CHIABSN = zeros(length(FAV), 1);
%     KSS =  zeros(length(FAV), 1);
%     for ii = 1:length(FAV)
%         FA = FAV(ii);
% 
%         [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
%         CHIABSN(ii) = CHIABS*N;
%         KSS(ii) = KS;
%     end
%     save(filename1, 'N', 'FAV', 'CHIABSN', 'PHIP', 'KSS'); 


    % Figure 1A
    N = 1e3;
    PHIPV = [0.2:0.1:0.9, 1-1e-4];
    CHIABSN = zeros(length(FAV), length(PHIPV));
    KSS =  zeros(length(FAV), length(PHIPV));

    for ii = 1:length(PHIPV)
        PHIP = PHIPV(ii)
        for jj = 1:length(FAV)
            FA = FAV(jj);
            [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
            CHIABSN(ii, jj) = CHIABS*N;
            KSS(ii, jj) = KS;
        end
    end
    save(filename1, 'N', 'FAV', 'CHIABSN', 'PHIPV', 'KSS'); 

% 
%     % Figure 1B
%     N = 1e0;
%     PHIPV = [0.2:0.1:0.9, 1-1e-4];
%     CHIABSN = zeros(length(FAV), length(PHIPV));
%     KSS =  zeros(length(FAV), length(PHIPV));
% 
%     for ii = 1:length(PHIPV)
%         PHIP = PHIPV(ii)
%         for jj = 1:length(FAV)
%             FA = FAV(jj);
%             [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
%             CHIABSN(ii, jj) = CHIABS*N;
%             KSS(ii, jj) = KS;
%         end
%     end
%     save(filename2, 'N', 'FAV', 'CHIABSN', 'PHIPV', 'KSS'); 

else
%     
%     % Figure 1A
%     load(filename1)
%     
%     figure;
%     plot(FAV, CHIABSN*PHIP, 'k-', 'linewidth', 2)
%     xlabel('f_A');ylabel('\phi_P\chi_{AB}^*N')
%     axis([0,1,0,80])
%     box on
%     set(gca,'fontsize',18)
% 
%     figure;
%     plot(FAV, KSS, 'k-', 'linewidth', 2)
%     xlabel('f_A');ylabel('q^*')
%     box on
%     set(gca,'fontsize',18)



    % Figure 1A
    load(filename1)
    
    figure;hold;
    for ii = 1:length(PHIPV)
        PHIP = PHIPV(ii);
        COL = (ii - 1) / (length(PHIPV) - 1);
        plot(FAV, CHIABSN(ii,:)*PHIP, '-', 'linewidth', 2, 'color', [COL 0 1-COL])
    end
    xlabel('f_A');ylabel('\phi_P\chi_{AB}^*vN')
    axis([0,1,0,80])
    box on
    set(gca,'fontsize',18)
    axis([0, 1, 0, 30])

    figure;hold
    for ii = 1:length(PHIPV)
        COL = (ii - 1) / (length(PHIPV) - 1);
        plot(FAV, KSS(ii,:)*sqrt(r2(N)), '-', 'linewidth', 2, 'color', [COL 0 1-COL])
    end
    xlabel('f_A');ylabel('q^*')
    box on
    set(gca,'fontsize',18)
    
    % Figure 1B
    load(filename2)
    
    figure;hold;
    for ii = 1:length(PHIPV)
        PHIP = PHIPV(ii);
        COL = (ii - 1) / (length(PHIPV) - 1);
        plot(FAV, CHIABSN(ii,:)*PHIP, '-', 'linewidth', 2, 'color', [COL 0 1-COL])
    end
    xlabel('f_A');ylabel('\phi_P\chi_{AB}^*vN')
    axis([0,1,0,80])
    box on
    set(gca,'fontsize',18)
    axis([0, 1, 0, 30])

    figure;hold
    for ii = 1:length(PHIPV)
        COL = (ii - 1) / (length(PHIPV) - 1);
        plot(FAV, KSS(ii,:)*sqrt(r2(N)), '-', 'linewidth', 2, 'color', [COL 0 1-COL])
    end
    xlabel('f_A');ylabel('q^*')
    box on
    set(gca,'fontsize',18)
end