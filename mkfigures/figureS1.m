clear;
addpath('../functions/')
addpath(genpath('../chainstats'))

figure1 = 0;
figure2 = 1;
figure3 = 0;

N = 1e0;
RM = sqrt(r2(N));
KV = logspace(-1, 2, 50)/RM;

if figure1
    PHIP = 1-1e-4;
    FAV = linspace(.1, .9, 101);
    CHIABSN = zeros(length(FAV), 1);
    for ii = 1:length(FAV)
        FA = FAV(ii);

        [CHIABS, KS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP);
        CHIABSN(ii) = CHIABS*N;

        [FA, CHIABS*N]
    end
    figure;plot(FAV, CHIABSN*PHIP, 'linewidth', 2)
end

if figure2
%     CHIV = linspace(0,2e3,51);
    CHIV = linspace(0,20,51);
    IEIG = 1;
    
    hfig = figure;
    loops = length(CHIV);
    F(loops) = struct('cdata',[],'colormap',[]);
    for j = 1:loops
        j
        
        CHIAB = CHIV(j);
        CHI = [0,CHIAB;CHIAB,0]/N;
%         CHI = [CHIAB,0;0,CHIAB]/N;
        
        plotphase_wsolvent_load(N, KV, CHI, IEIG)
        F(j) = getframe(hfig);
    end
    movie2avi(F,'test.avi')
end

if figure3
    FA = 0.5;PHIP = 0.5;
%     [CHIABS, KS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP);
    
    figure;hold
    cnt = 1;
    for CHIAB = [0, 10]
        COL = (cnt -1) / (2-1)
        CHIBA = CHIAB;
        CHI = [0, CHIAB; CHIBA, 0]/N;

        [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);

        plot(KV, 1./EIG(:,1)/N, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
        plot(KV, 1./EIG(:,2)/N, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
        cnt = cnt + 1;
    end

    box on
    xlabel('2l_Pq');ylabel('<\psi(-q)\psi(q)>')
    xlim([min(KV), max(KV)])
    set(gca,'fontsize',18)
    set(gca,'xscale','log')
    set(gca,'yscale','log')

end

% 
% CHIABSN = []
% for N = logspace(-1, 4, 20)
%     [CHIABS, KS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP);
%     CHIABSN = [CHIABSN, CHIABS*N];
%     
%     [N, CHIABSN(end)]
% end
% figure;semilogx(logspace(-1, 3, 20), CHIABSN*PHIP, 'linewidth', 2)

% % [CHIABS,KSS,EIGVS] = spinodal_wsolvent(N,FA,KV,PHIP);
% figure;hold
% CHIABS = spinodal_wsolvent(N,FA,KV,PHIP);


% CHI=0;
% val=gamma2(N,FA,kv,CHI);
% figure;
% plot(kv*N, 1./val)
% set(gca,'xscale','log')
% set(gca,'yscale','log')

% 
% s2invAA = [];
% s2AA = [];
% for k = kv
%     s2inv = s2inverse(N,FA,k);
%     s2 = s2gc(N,FA,k);
%     s2invAA = [s2invAA, s2inv(1, 1)];
%     s2AA = [s2AA, s2(1, 1)];
% end
% 
% figure;loglog(kv*N, s2invAA)
% figure;loglog(kv*N, s2AA)