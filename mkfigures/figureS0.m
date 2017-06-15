clear;

PHIP = 0.5;
for N = [1e3, 1e0];
%     for FA = [0.5, 0.8];
    for FA = [0.2]

    % calculate spinodal CHIABS
    [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);

    RM = sqrt(r2(N));
    KV = logspace(-1, 4, 1000)/RM;

    CHIABV = fliplr([0,0.2,0.4,0.6,0.8]);
    figure;hold
    for ii = 1:length(CHIABV)
        COL = 1-(ii-1)/(length(CHIABV)-1)

        CHIAB = CHIABV(ii)*CHIABS;
        CHIBA = CHIAB;
        CHI = [0, CHIAB; CHIBA, 0];

        [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);
        plot(KV*RM, 1./EIG(:,1)/N, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
    end
    
    for ii = 1:length(CHIABV)
        COL = 1-(ii-1)/(length(CHIABV)-1)

        CHIAB = CHIABV(ii)*CHIABS;
        CHIBA = CHIAB;
        CHI = [0, CHIAB; CHIBA, 0];

        [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);
        plot(KV*RM, 1./EIG(:,2)/N, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
    end
    
    % plot again
    CHIABV = [0,0.2,0.4,0.6,0.8];
    for ii = 1:length(CHIABV)
        COL = (ii-1)/(length(CHIABV)-1)

        CHIAB = CHIABV(ii)*CHIABS;
        CHIBA = CHIAB;
        CHI = [0, CHIAB; CHIBA, 0];

        [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);
        plot(KV*RM, 1./EIG(:,1)/N, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
    end
    
    for ii = 1:length(CHIABV)
        COL = (ii-1)/(length(CHIABV)-1)

        CHIAB = CHIABV(ii)*CHIABS;
        CHIBA = CHIAB;
        CHI = [0, CHIAB; CHIBA, 0];

        [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);
        plot(KV*RM, 1./EIG(:,2)/N, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
    end

    box on
    xlabel('Rq');
    ylabel('$\langle {\psi}^{\bf{\top}}(\vec{q}) {\psi}(-\vec{q}) \rangle/Nv$',...
        'Interpreter', 'Latex')
    xlim([min(KV), max(KV)]*RM)
    set(gca,'fontsize',18)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
%     title(strcat('N=', sprintf('%d', N), ', f_A=', sprintf('%.2f', FA)))

    h_leg=legend('\chi_{AB}=0.8\chi_{AB}^{*}','\chi_{AB}=0.6\chi_{AB}^{*}',  ...
        '\chi_{AB}=0.4\chi_{AB}^{*}',...
        '\chi_{AB}=0.2\chi_{AB}^{*}', '\chi_{AB}=0','location', 'northeast');
    set(h_leg,'box','off')
    axis([1e-1,1e4,1e-6,1e2])
    set(gca,'ytick',[1e-6,1e-4,1e-2,1e0,1e2])
    set(gca,'xtick',[1e-1,1e0,1e1,1e2,1e3,1e4])
    if N==1e3
        plotlog(5e2,1e-4,-1,2)
    else
        plotlog(5e2,2e-3,-1,2)
    end
    end
end