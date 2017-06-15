% Figure 3: Spinodals with changing solvent concentration
clear;

CALCON = 0;
filename1 = '../data/figureS3A';
filename2 = '../data/figureS3B';

if CALCON
    
    % Figure 3A
    PHIPV = logspace(-1,-0.00001,101);
    NV = logspace(0,3,7);
    
    FA = 0.5;
    CHIABSN = zeros(length(NV), length(PHIPV));
    KSS =  zeros(length(NV), length(PHIPV));
    
    for ii = 1:length(NV)
        N = NV(ii)

        for jj = 1:length(PHIPV)
            PHIP = PHIPV(jj);

            [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
            CHIABSN(ii, jj) = CHIABS*N;
            KSS(ii, jj) = KS;
        end
    end
    save(filename1, 'NV', 'PHIPV', 'CHIABSN', 'KSS'); 
    
    % Figure 3B and C
    FA = 0.2;
    CHIABSN = zeros(length(NV), length(PHIPV));
    KSS =  zeros(length(NV), length(PHIPV));
    for ii = 1:length(NV)
        N = NV(ii)

        for jj = 1:length(PHIPV)
            PHIP = PHIPV(jj);

            [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
            CHIABSN(ii, jj) = CHIABS*N;
            KSS(ii, jj) = KS;
        end
    end
    save(filename2, 'NV', 'PHIPV', 'CHIABSN', 'KSS'); 

else
    close all
    
%     % Figure 3A
%     load(filename1)
%     
%     figure(1);hold
%     for ii = 1:length(NV)
%         COL = (ii-1)/(length(NV)-1);
%         figure(1);plot(PHIPV, CHIABSN(ii,:), 'color', [COL 0 1-COL], 'linewidth', 2)
%     end
%     
%     figure(1);
%     box on
%     xlabel('\phi_P');ylabel('\chi_{AB}^*N')
%     set(gca,'xscale','log');set(gca,'yscale','log')
%     set(gca,'fontsize',18)
%     % x = logspace(-.3,0,3);y = 28*power(x, -1);
%     % plot(x, y, 'k--', 'linewidth', 2)
%     set(gca,'xtick',[0.1, 0.3, 0.6, 1.0])
%     set(gca,'xticklabels', {'0.1', '0.3', '0.6', '1.0'})
%     ylim([5,1000])

    
    
    % Figure 3A
    load(filename1)
    
    figure(1);hold
    figure(2);hold
    for ii = 1:length(NV)
        N = NV(ii)
        COL = (ii-1)/(length(NV)-1);
        figure(1);plot(PHIPV, CHIABSN(ii,:), 'color', [COL 0 1-COL], 'linewidth', 2)
        figure(2);plot(PHIPV, KSS(ii,:)*sqrt(r2(N)), 'color', [COL 0 1-COL], 'linewidth', 2)
    end

    figure(1);
    box on
    xlabel('\phi_P');ylabel('\chi_{AB}^*N')
    set(gca,'xscale','log');set(gca,'yscale','log')
    set(gca,'fontsize',18)
    % x = logspace(-.3,0,3);y = 28*power(x, -1);
    % plot(x, y, 'k--', 'linewidth', 2)
    set(gca,'xtick',[0.1, 0.3, 0.6, 1.0])
    set(gca,'xticklabels', {'0.1', '0.3', '0.6', '1.0'})
    ylim([5,1000])

    figure(2);
    box on
    xlabel('\phi_P');ylabel('Rq_1^*(\chi_{AB}^*)')
    set(gca,'xscale','log');
    set(gca,'fontsize',18)
    set(gca,'xtick',[0.1, 0.3, 0.6, 1.0])
    set(gca,'xticklabels', {'0.1', '0.3', '0.6', '1.0'})

    set(gca,'ytick',[5.0, 5.4, 5.8, 6.2, 6.6])
    set(gca,'yticklabels', {'5.0', '5.4', '5.8', '6.2', '6.6'})
    
    
    
    
    
    % Figure 3B
    load(filename2)
    
    figure(3);hold
    figure(4);hold
    for ii = 1:length(NV)
        N = NV(ii)
        COL = (ii-1)/(length(NV)-1);
        figure(3);plot(PHIPV, CHIABSN(ii,:), 'color', [COL 0 1-COL], 'linewidth', 2)
        figure(4);plot(PHIPV, KSS(ii,:)*sqrt(r2(N)), 'color', [COL 0 1-COL], 'linewidth', 2)
    end

    figure(3);
    box on
    xlabel('\phi_P');ylabel('\chi_{AB}^*N')
    set(gca,'xscale','log');set(gca,'yscale','log')
    set(gca,'fontsize',18)
    % x = logspace(-.3,0,3);y = 28*power(x, -1);
    % plot(x, y, 'k--', 'linewidth', 2)
    set(gca,'xtick',[0.1, 0.3, 0.6, 1.0])
    set(gca,'xticklabels', {'0.1', '0.3', '0.6', '1.0'})
    ylim([5,1000])

    figure(4);
    box on
    xlabel('\phi_P');ylabel('Rq_1^*(\chi_{AB}^*)')
    set(gca,'xscale','log');
    set(gca,'fontsize',18)
    set(gca,'xtick',[0.1, 0.3, 0.6, 1.0])
    set(gca,'xticklabels', {'0.1', '0.3', '0.6', '1.0'})

    set(gca,'ytick',[5.0, 5.4, 5.8, 6.2, 6.6])
    set(gca,'yticklabels', {'5.0', '5.4', '5.8', '6.2', '6.6'})
end