clear;
addpath('../functions/')
addpath(genpath('../chainstats'))

figure1 = 0;
figure2 = 1;
figure3 = 0;
figure4 = 0;
figure5 = 0;


N = 10;
FAV = linspace(0, 0.01, 10);
CHIABS = zeros(10, 1);
PHIP = 0.5;

cnt = 1;
for FA = FAV
    CHIABS(cnt) = spinodal_wsolvent(N, FA, PHIP);
    [cnt, CHIABS(cnt)]
    cnt = cnt + 1;
end


if figure1
    N = 1e4;
    FAV = linspace(.01, .99, 21);
    
    figure;hold
    cnt = 1;
    LOAD = 1;
    FAIND = 1;
    
    NV = [1e-1,1e0,1e1,1e2,1e3];
    for N = NV
        COL = (cnt-1)/(length(NV)-1);
        
        [FAV, PHIPV, CHISNPHIP, KSS, EIGVSS] = plotspin_wsolvent(N, LOAD);
        plot(CHISNPHIP(FAIND,:)./PHIPV, PHIPV, 'linewidth', 2, 'color', [COL 0 1-COL])

        cnt = cnt+1;
    end
    xlabel('\chiN');ylabel('\phi_{ODT}')
    box on
    set(gca,'fontsize',18)
    set(gca,'xscale','log');set(gca,'yscale','log')
    axis([5,500,0.1,1])
    legend('N=1','N=10','N=100','N=1000')

end

if figure2
    LOAD = 1;
    
    NV = fliplr(logspace(-1, 3, 9));
%     for N = NV
%         CHIV = linspace(0,20,51)/N;
%         IEIG = 1;
%         loops = length(CHIV);
% 
%         hfig = figure;
%         F(loops) = struct('cdata',[],'colormap',[]);
%         for j = 1:loops
%             j
% 
%             CHIAB = CHIV(j);
%             CHI = [0,CHIAB;CHIAB,0];
% 
%             [FAV, PHIPV, EIGV, EIG, KSV] = plotphase_wsolvent(N, CHI, IEIG, LOAD);
%             F(j) = getframe(hfig);
%         end
%         movie2avi(F,'test.avi')
%     end
    
    NV = [1e0,1e1,1e2,1e3];
    cnt = 1;
%     figure;hold
    figure(3);
    for N = NV
        N
        
        CHIV = linspace(0,20,51)/N;
        IEIG = 1;
        KSS = zeros(length(CHIV),1);
        loops = length(CHIV);

        for j = 1:loops
            CHIAB = CHIV(j);
            CHI = [0,CHIAB;CHIAB,0];

            [FAV, PHIPV, EIGV, EIG, KSV] = plotphase_wsolvent(N, CHI, IEIG, LOAD);
            KSS(j) = KSV(21,26);
        end
        
        COL = (cnt-1) / (length(NV)-1);
        plot(CHIV*N, KSS, '--','linewidth', 2, 'color', [COL 0 1-COL])
        cnt = cnt + 1;
    end
    
    box on
    set(gca,'fontsize',18)
    xlabel('\chiN');ylabel('K*')
    legend('N=1','N=10','N=100','N=1000')
end

if figure3
    N = 1e0;
    RM = sqrt(r2(N));
    KV = logspace(-1, 2, 100)/RM;

    FA = 0.5;PHIP = 0.9999;
    [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP)
    
    figure;hold
    cnt = 1;NCHIAB=3;
    for CHIAB = linspace(0,10,NCHIAB)
        COL = (cnt - 1) / (NCHIAB - 1)
        CHIBA = CHIAB;
        CHI = [0, CHIAB; CHIBA, 0]/N/PHIP;

        [EIG,EIGV]=gamma2_solvent(N,FA,KV,CHI,PHIP);

        plot(KV*RM, 1./EIG(:,1)/N, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
        plot(KV*RM, 1./EIG(:,2)/N, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
        cnt = cnt + 1;
    end

    box on
    xlabel('Rq');ylabel('<\psi(-q)\psi(q)>')
    xlim([min(KV), max(KV)]*RM)
    set(gca,'fontsize',18)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end

if figure4
    LOAD = 0;
    NV = fliplr(logspace(-1, 3, 9));
    NV = [0.1]
    for N = NV
        [FAV, PHIPV, CHISNPHIP, KSS] = plotspin_wsolvent(N, LOAD);
    end
end

if figure5
    CHIABSN = []
    NV = fliplr(logspace(-1, 3, 9));
    for N = NV
        [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);
        CHIABSN = [CHIABSN, CHIABS*N*PHIP];

        [N, CHIABSN(end)]
    end
    figure;semilogx(NV, CHIABSN*PHIP, 'linewidth', 2)
    
end

% CHI=0;
% val=gamma2(N,FA,kv,CHI);
% figure;
% plot(kv*N, 1./val)
% set(gca,'xscale','log')
% set(gca,'yscale','log')

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