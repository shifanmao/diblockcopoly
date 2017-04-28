clear;
addpath('../functions/')
addpath(genpath('../chainstats'))

figure1 = 0;
figure2 = 0;
figure3 = 0;
figure4 = 1;

N = 1e4;
RM = sqrt(r2(N));
KV = logspace(-1, 2, 100)/RM;

if figure1
    figure;hold
    
    PHIP = 1-1e-4;
    FAV = linspace(.1, .9, 51);
    CHIABSN = zeros(length(FAV), 1);
    for ii = 1:length(FAV)
        FA = FAV(ii);

        [CHIABS, KS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP);
        CHIABSN(ii) = CHIABS*N;

        [FA, CHIABS*N]
    end
    
    plot(FAV, CHIABSN*PHIP, 'k','linewidth', 2)
    
    PHIP = 0.5;
    FAV = linspace(.1, .9, 51);
    CHIABSN = zeros(length(FAV), 1);
    for ii = 1:length(FAV)
        FA = FAV(ii);

        [CHIABS, KS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP);
        CHIABSN(ii) = CHIABS*N;

        [FA, CHIABS*N]
    end
    
    plot(FAV, CHIABSN*PHIP, 'kx','linewidth', 2)
    
    xlabel('f_A');ylabel('\chi_{spin}N\phi_P')
    xlim([FAV(1),FAV(end)])
    box on
    set(gca,'fontsize',18)
end

if figure2
   map = ...
        [0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0];
    
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
        colormap(map)
        F(j) = getframe(hfig);
    end
    movie2avi(F,'test.avi')
end

if figure3
    FA = 0.5;PHIP = 0.9999;
    [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, KV, PHIP)
    
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
    
    FAV = linspace(.1, .9, 21);
    PHIPV = linspace(.001, .999, 11);
    PHIPV = .9999;
    CHISNPHIP = zeros(length(FAV), length(PHIPV));
    KSS = zeros(length(FAV), length(PHIPV));

    % filename='data/newgamdata';
    % outfile = fopen(filename, 'wt');
    for jj = 1:length(PHIPV)
        PHIP = PHIPV(jj)
        for ii = 1:length(FAV)
            FA = FAV(ii);

%             CHIS(ii,jj) = spinodal_wsolvent(N, FA, KV, PHIP);
            [CHIABS, KS, ~, ~] = spinodal_wsolvent(N, FA, KV, PHIP);
            CHISNPHIP(ii,jj)=CHIABS*N*PHIP;
            KSS(ii,jj)=KS;
        end
    end

%     figure;plotcontour(FAV, 1-PHIPV, CHISNPHIP)
    figure;plotsurf(FAV, 1-PHIPV, CHISNPHIP)
    figure;plotsurf(FAV, 1-PHIPV, KSS)
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