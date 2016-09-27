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
for EPS = 1
    figure;hold;set(gca,'fontsize',20)
    
    % load simulation paramters
    simparams = load('test/simparams');
    N =EPS*Gtimes2;
    C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);C = sqrt(C2);
    
    RUNNUM = 2;LBOX = 20;
    PREF = plotsimfig(EPS,Gtimes2,RUNNUM,simparams,LBOX,DEL)
    
    RUNNUM = 3;LBOX = 16;
    PREF = plotsimfig(EPS,Gtimes2,RUNNUM,simparams,LBOX,DEL)
    
    PLOTDENSITY = 0;PLOTRG = 1;NCHI = 50;
    densityRG(N,C,0.5,PLOTDENSITY,PLOTRG,NCHI);
    legend('MC Sim. L_{box}=20','MC Sim. L_{box}=16','Mean-field Theory','F-H Theory')
    
    xlim([0,20]);
    title(sprintf('N=%.2f, C^2=%.2f',N,C^2))

%     % imagename = sprintf('mkfigures/SINV_EPS%.2fG%dLBOX%d.png',EPS,Gtimes2,LBOX);
%     % saveas(gca,imagename,'png')
%     imagename = sprintf('mkfigures/SINV_N%.1fG%d_PREF1.eps',N,Gtimes2);
%     saveas(gca,imagename,'eps')
end

cd mkfigures/
end

function PREF = plotsimfig(EPS,Gtimes2,RUNNUM,simparams,LBOX,DEL)

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

% % correction according to Matsen find RPA density correlations

% assume as RPA density fluctuations
K = 0:(LBOX/DEL);
k = [];
for IX = K
    for IY = K
        for IZ = K
            if norm([IX,IY,IZ])>0
                k = [k,2*pi/LBOX*(2*LP)*norm([IX,IY,IZ])];
            end
        end
    end
end
SRPA = (1/V)./gamma2(N,FA,k,0);
A = (2*pi)^3/LBOX^3*sum(SRPA)

% % alternatively load density fluctuations from simulation
% SAVEFILENAME = sprintf('SDIB_SIM_CHIG%.3fEPS%.2fG%dFA%.2f',0.08,EPS,Gtimes2,FA);
% SIMdata = load(strcat('~shifan/Documents/research/data-diblock-09-17-16/',filename,...
%     '/mkfigures/savedata/',SAVEFILENAME));
% KSIM = SIMdata(:,1);
% SSIM = (1/V)*EPS*mean(SIMdata(:,2:end),2);
% SRPA = 4*pi*KSIM(1:end-1).^2.*(KSIM(2:end)-KSIM(1:end-1)).*SSIM(1:end-1);
% A = sum(SRPA)

% calculate prefactor
PREF = 1-A*V^2/((2*pi)^3*FA*(1-FA));

% % make a plot of integrated density fluctuations
% kmin = 2*pi/LBOX*(2*LP)*1;
% kmax = 2*pi/LBOX*(2*LP)*20;
% k = logspace(log10(kmin),log10(kmax),50);
% SRPA_PLOT = (1/V)./gamma2(N,FA,k,0);

% figure;hold
% plot(k,SRPA_PLOT)
% plot(KSIM*(2*LP),SSIM)

% % test behavior of SRPA
% N = 1e3;
% k = logspace(-1,1.5,100)';
% R = sqrt(r2(N));
% SRPA = (1/V)./gamma2(N,FA,k,0);
% Y = SRPA*V/N - 12*FA*(1-FA)./power(k*R,2);
% figure;loglog(k,Y)

% % correction according to Olvera de la Cruz
% PREF = 1-6*1.221/pi^2*0.1/L0^2*pi;

sinv = load(strcat('~shifan/Documents/research/data-diblock-09-17-16/',filename,...
    '/mkfigures/savedata/SDIB_SIM_SINV'));

if LBOX == 20
    plot(sinv(:,1)*PREF,sinv(:,2)*Gtimes2,'b^','markersize',8)
elseif LBOX == 16
    plot(sinv(:,1)*PREF,sinv(:,2)*Gtimes2,'b<','markersize',8)
end
% errorbar(sinv(:,1),sinv(:,2)*Gtimes2,sinv(:,3)*Gtimes2,'bo','markersize',8)
end