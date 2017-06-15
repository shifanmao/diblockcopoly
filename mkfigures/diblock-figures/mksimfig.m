function mksimfig
cd ../

clear;
%close all;
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

DEL = 1;
Gtimes2 = 32;
% for EPS = [1.0,0.5,0.1,0.05]
for EPS = [0.05]
figure('Position', [100, 100, 1200, 900]);
hold;set(gca,'fontsize',50);
    
    % load simulation paramters
%     simparams = load('test/simparams');
    simparams = load(strcat('~shifan/Documents/research/data-diblock-09-17-16/',...
        'utilitycode/utility/sdata/simparams'));
    N =EPS*Gtimes2;

    RUNNUM = 2;LBOX = 20;
    if (RUNNUM == 2 || RUNNUM == 3)
        C = sqrt(1525.9);
    end
    [PREF1,KRPA1,SRPA1,CHISIM1,SINVSIM1] = plotsimfig(EPS,Gtimes2,RUNNUM,simparams,LBOX,DEL);
    
    RUNNUM = 3;LBOX = 16;
    if (RUNNUM == 2 || RUNNUM == 3)
        C = sqrt(1525.9);
    end
    [PREF2,KRPA2,SRPA2,CHISIM2,SINVSIM2] = plotsimfig(EPS,Gtimes2,RUNNUM,simparams,LBOX,DEL);
    
    PLOTDENSITY = 0;PLOTRG = 1;NCHI = 50;
    densityRG(N,C,0.5,0,PLOTDENSITY,PLOTRG,NCHI);
    
    % plot simulation results
    SIMEPS = 32*[0.05,0.1,0.5,1.0];
    SIMODT1 = [16.361,15.958,15.914,16.094];
    SIMODT2 = [16.672,16.533,16.311,16.374];
    SIMODT = (SIMODT1+SIMODT2)/2;
    IND = find(SIMEPS==N);
    plot([SIMODT1(IND),SIMODT1(IND)]*PREF1,[0,20],'c-','linewidth',2)
    plot([SIMODT2(IND),SIMODT2(IND)]*PREF1,[0,20],'m-','linewidth',2)
        
    xlim([0,25]);
    legend('L_{box}=20','L_{box}=16')
%     legend('MC Sim. L_{box}=20','MC Sim. L_{box}=16','Mean-field Theory','F-H Theory')
%    title(sprintf('N=%.2f, C^2=%.2f',N,C^2))
    imagename = sprintf('mkfigures/SINV_N%.1fG%d.eps',N,Gtimes2);
    saveas(gca,imagename,'epsc')
end

cd mkfigures/
end

function [PREF,KRPA,SRPA,CHISIM,SINVSIM] = plotsimfig(EPS,Gtimes2,RUNNUM,simparams,LBOX,DEL)

% parse EPS
if EPS==1
    filename = sprintf('EPS1-N%d-RUN0%d',Gtimes2,RUNNUM);
elseif EPS == 0.5
    filename = sprintf('EPSp5-N%d-RUN0%d',Gtimes2,RUNNUM);
elseif EPS == 0.1
    filename = sprintf('EPSp1-N%d-RUN0%d',Gtimes2,RUNNUM);
elseif EPS == 0.05
    filename = sprintf('EPSp05-N%d-RUN0%d',Gtimes2,RUNNUM);
elseif EPS == 0.01
    filename = sprintf('EPSp01-N%d-RUN0%d',Gtimes2,RUNNUM);
end

% load simulation parameters
N =EPS*Gtimes2;
FA = 0.5;
V = 0.1;
LP = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,3)/2;
C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);C = sqrt(C2);

% % % correction according to Matsen find RPA density correlations
% % % assume as RPA density fluctuations
% K = 0:(LBOX/DEL);
% KRPA = [];
% for IX = K
%     for IY = K
%         for IZ = K
%             if norm([IX,IY,IZ])>0
%                 KRPA = [KRPA,2*pi/LBOX*(2*LP)*norm([IX,IY,IZ])];
%             end
%         end
%     end
% end
% SRPA = (1/V)./(EPS*gamma2(N,FA,KRPA,0));
% A = (2*pi)^3/LBOX^3*sum(SRPA);
% 
% % calculate prefactor
% PREF = 1-A*V^2/((2*pi)^3*FA*(1-FA))
KRPA = 1;
SRPA = 1;
PREF = 1;

% % % alternatively load density fluctuations from simulation
% % SAVEFILENAME = sprintf('SDIB_SIM_CHIG%.3fEPS%.2fG%dFA%.2f',0.08,EPS,Gtimes2,FA);
% % SIMdata = load(strcat('~shifan/Documents/research/data-diblock-09-17-16/',filename,...
% %     '/mkfigures/savedata/',SAVEFILENAME));
% % KSIM = SIMdata(:,1);
% % % SSIM = (1/V)*EPS*mean(SIMdata(:,2:end),2);
% % SSIM = (1/V)*mean(SIMdata(:,2:end),2);
% % SRPA = 4*pi*KSIM(1:end-1).^2.*(KSIM(2:end)-KSIM(1:end-1)).*SSIM(1:end-1);
% % A = sum(SRPA);
% % 


% % make a plot of integrated density fluctuations
% kmin = 2*pi/LBOX*(2*LP)*1;
% kmax = 2*pi/LBOX*(2*LP)*20;
% k = logspace(log10(kmin),log10(kmax),50);
% SRPA_PLOT = (1/V)./(EPS*gamma2(N,FA,k,0));
% 
% figure;hold
% plot(k,SRPA_PLOT)
% plot(KSIM*(2*LP),SSIM)
% set(gca,'xscale','log');set(gca,'yscale','log')

% % test behavior of SRPA
% NV = logspace(4,5,10);
% Y1 = [];
% Y2 = [];
% for N = NV
% % k = logspace(-1,1.5,100)';
% R = sqrt(r2(N));V = 1;FA=0.5;
% k = 1e9;
% SRPA = (1/V)./(gamma2(N,FA,k,0));
% % Y = SRPA*V/N - 12*FA*(1-FA)./power(k*R,2);
% Y1 = [Y1,SRPA*V/N];
% Y2 = [Y2,12*FA*(1-FA)./power(k*R,2)];
% % Y = [Y,SRPA*V];
% end
% figure;hold
% plot(NV,Y1,'-');
% plot(NV,Y2,'x');
% legend('y1','y2')
% set(gca,'xscale','log');set(gca,'yscale','log')

% % correction according to Olvera de la Cruz
% PREF = 1-6*1.221/pi^2*0.1/L0^2*pi;

sinv = load(strcat('~shifan/Documents/research/data-diblock-09-17-16/',filename,...
    '/mkfigures/savedata/SDIB_SIM_SINV'));
CHISIM = sinv(:,1);
SINVSIM = sinv(:,2);

if LBOX == 20
    whos CHISIM
    plot(CHISIM*PREF,SINVSIM*Gtimes2,'k^','markersize',18)
elseif LBOX == 16
    plot(CHISIM*PREF,SINVSIM*Gtimes2,'k<','markersize',18)
end
% errorbar(sinv(:,1),sinv(:,2)*Gtimes2,sinv(:,3)*Gtimes2,'bo','markersize',8)
end