% Figure 2: Triangular phases at 8different 
clear;

CALCON = 0;
filename1 = '../data/figureS2ATriMesh';
filename2 = '../data/figureS2BTriMesh';
filename3 = '../data/figureS2CTriMesh';
filename4 = '../data/figureS2DTriMesh';
filenames = {filename1, filename2, filename3, filename4};

if CALCON
    % Figure 2A-D
    cnt = 1;
    for N = [1e0, 1e3]
        for CHI = [0, 0.8]
            filename = filenames{cnt};
            
            [FAVV, PHIPV, EIGV, EIG, KSV] = calcphase_wsolvent(N, CHI, 1, 0);
            save(filename, 'FAVV', 'PHIPV', 'EIGV', 'EIG', 'KSV', 'N');
            cnt = cnt+1;
        end
    end
    
else
    for ii = 1:4
        filename = filenames{ii};
        load(filename);
        
        figure;plotphase_wsolvent(FAVV, PHIPV, EIGV, EIG, KSV, N);
        if ii == 1 || ii == 3
            title('\chi_{AB}=0');set(get(gca,'title'), 'position', [0.2,1])
        else
            title('\chi_{AB}=0.8\chi_{AB}^*');set(get(gca,'title'), 'position', [0.2,1])
        end
    end
    
    PHIP = 0.5;
    for N = [1e3, 1e0];
        for FA = [0.5, 0.8];

        % calculate spinodal CHIABS
        [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);

        RM = sqrt(r2(N));
        KV = logspace(-1, 2, 100)/RM;

        CHIABV = [0,0.2,0.4,0.6,0.8];
        figure;hold
        for ii = 1:length(CHIABV)
            COL = (ii-1)/(length(CHIABV)-1)

            CHIAB = CHIABV(ii)*CHIABS;
            CHIBA = CHIAB;
            CHI = [0, CHIAB; CHIBA, 0];

            [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);

            plot(KV*RM, 1./EIG(:,1)/N, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
            plot(KV*RM, 1./EIG(:,2)/N, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
        end

        box on
        xlabel('Rq');
%         ylabel('<\psi(-q)\psi(q)>')
        ylabel('$\langle \bf{\psi}^{\bf{\top}}(\vec{q}) \bf{\psi}(-\vec{q}) \rangle/Nv$',...
            'Interpreter', 'Latex')
        xlim([min(KV), max(KV)]*RM)
        set(gca,'fontsize',18)
        set(gca,'xscale','log')
        set(gca,'yscale','log')
%         title(strcat('N=', sprintf('%d', N), ', f_A=', sprintf('%.2f', FA)))

        if N==1e3
            legend('\chi_{AB}=0', '\chi_{AB}=0.2\chi_{AB}^{*}', '\chi_{AB}=0.4\chi_{AB}^{*}',...
           '\chi_{AB}=0.6\chi_{AB}^{*}', '\chi_{AB}=0.8\chi_{AB}^{*}', 'location', 'northwest')
        else
            legend('\chi_{AB}=0', '\chi_{AB}=0.2\chi_{AB}^{*}', '\chi_{AB}=0.4\chi_{AB}^{*}',...
           '\chi_{AB}=0.6\chi_{AB}^{*}', '\chi_{AB}=0.8\chi_{AB}^{*}', 'location', 'southeast')
        end
        end
    end
end