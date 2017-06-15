% Figure 1: Spinodals of diblock copolymer solutions
clear;

% Figure 1: Spinodals of diblock copolymer solutions
clear;

CALCON = 0;
filename1 = '../data/figureS1A';
filename2 = '../data/figureS1B';
filenames = {filename1, filename2};

if CALCON
    FAV = linspace(.1, .9, 101);
    PHIPV = [0.1:0.1:1];
    
    NV = [1e3, 1e0];
    for ff = 1:2
        N = NV(ff);
        filename = filenames{ff};
        
        CHIABSN = zeros(length(PHIPV), length(FAV));
        KSS =  zeros(length(PHIPV), length(FAV));

        for ii = 1:length(PHIPV)
            PHIP = PHIPV(ii);
            for jj = 1:length(FAV)
                FA = FAV(jj);
                [PHIP, FA]
                [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
                CHIABSN(ii, jj) = CHIABS*N;
                KSS(ii, jj) = KS;

%                 [KS, ~, ~] = eigs_wsolvent(N, FA, PHIP, [0,0;0,0], 1);
%                 KSS(ii, jj) = KS;
            end
        end
        save(filename, 'N', 'FAV', 'CHIABSN', 'PHIPV', 'KSS'); 
    end
else
    for ff = 1:2
        filename = filenames{ff};
        load(filename)
        FAV = linspace(.1, .9, 101);
        PHIPV = [0.1:0.1:1];


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
        axis([.1, .9, 5, 20])

        figure;hold
        for ii = 1:length(PHIPV)
            COL = (ii - 1) / (length(PHIPV) - 1);
            plot(FAV, KSS(ii,:)*sqrt(r2(N)), '-', 'linewidth', 2, 'color', [COL 0 1-COL])
        end
        xlabel('f_A');ylabel('Rq^*(\chi_{AB}^*)')
        box on
        if ff==1
            axis([.1, .9, 4.5,5.5])
        else
            axis([.1, .9, 5.5, 7])
        end
        set(gca,'fontsize',18)
    end
end