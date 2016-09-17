cd ../

clear;close all;
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

simparams = load('~shifan/Desktop/data-diblock-08-31-06/utilitycode/utility/sdata/simparams');

EPS = 1;
Gtimes2 = 32;
N =EPS*Gtimes2;

C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);
C = sqrt(C2);

sinv = load('~shifan/Desktop/data-diblock-08-31-06/EPS1-N32-RUN01/mkfigures/savedata/SDIB_SIM_SINV');

PLOTRG = 1;NCHI = 100;
densityRG(N,C,0.5,PLOTRG,NCHI);
plot(sinv(:,1),sinv(:,2)*Gtimes2,'bo','markersize',8)
xlim([0,30])

imagename = sprintf('mkfigures/SINV_EPS%.2fG%d.png',EPS,Gtimes2);
saveas(gca,imagename,'png')


%%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%

EPS = 0.5;
Gtimes2 = 32;
N =EPS*Gtimes2;

C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);
C = sqrt(C2);

sinv = load('~shifan/Desktop/data-diblock-08-31-06/EPSp5-N32-RUN01/mkfigures/savedata/SDIB_SIM_SINV');

PLOTRG = 1;NCHI = 100;
densityRG(N,C,0.5,PLOTRG,NCHI);
plot(sinv(:,1),sinv(:,2)*Gtimes2,'bo','markersize',8)
xlim([0,30])

imagename = sprintf('mkfigures/SINV_EPS%.2fG%d.png',EPS,Gtimes2);
saveas(gca,imagename,'png')


%%%%%%%%%%%%%%% 3 %%%%%%%%%%%%%%%

EPS = 0.1;
Gtimes2 = 32;
N =EPS*Gtimes2;

C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);
C = sqrt(C2);

sinv = load('~shifan/Desktop/data-diblock-08-31-06/EPSp1-N32-RUN01/mkfigures/savedata/SDIB_SIM_SINV');

PLOTRG = 1;NCHI = 100;
densityRG(N,C,0.5,PLOTRG,NCHI);
plot(sinv(:,1),sinv(:,2)*Gtimes2,'bo','markersize',8)
xlim([0,30])

imagename = sprintf('mkfigures/SINV_EPS%.2fG%d.png',EPS,Gtimes2);
saveas(gca,imagename,'png')



%%%%%%%%%%%%%%% 4 %%%%%%%%%%%%%%%

EPS = 0.05;
Gtimes2 = 32;
N =EPS*Gtimes2;

C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);
C = sqrt(C2);

sinv = load('~shifan/Desktop/data-diblock-08-31-06/EPSp05-N32-RUN01/mkfigures/savedata/SDIB_SIM_SINV');

PLOTRG = 1;NCHI = 100;
densityRG(N,C,0.5,PLOTRG,NCHI);
plot(sinv(:,1),sinv(:,2)*Gtimes2,'bo','markersize',8)
xlim([0,30])

imagename = sprintf('mkfigures/SINV_EPS%.2fG%d.png',EPS,Gtimes2);
saveas(gca,imagename,'png')



%%%%%%%%%%%%%%% 5 %%%%%%%%%%%%%%%

EPS = 0.5;
Gtimes2 = 64;
N =EPS*Gtimes2;

C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);
C = sqrt(C2);

sinv = load('~shifan/Desktop/data-diblock-08-31-06/EPSp5-N64-RUN01/mkfigures/savedata/SDIB_SIM_SINV');

PLOTRG = 1;NCHI = 100;
densityRG(N,C,0.5,PLOTRG,NCHI);
plot(sinv(:,1),sinv(:,2)*Gtimes2,'bo','markersize',8)
xlim([0,30])

imagename = sprintf('mkfigures/SINV_EPS%.2fG%d.png',EPS,Gtimes2);
saveas(gca,imagename,'png')


%%%%%%%%%%%%%%% 6 %%%%%%%%%%%%%%%

EPS = 0.05;
Gtimes2 = 64;
N =EPS*Gtimes2;

C2 = simparams(simparams(:,1)==Gtimes2 & simparams(:,2) == EPS,4);
C = sqrt(C2);

sinv = load('~shifan/Desktop/data-diblock-08-31-06/EPSp05-N64-RUN01/mkfigures/savedata/SDIB_SIM_SINV');

PLOTRG = 1;NCHI = 100;
densityRG(N,C,0.5,PLOTRG,NCHI);
plot(sinv(:,1),sinv(:,2)*Gtimes2,'bo','markersize',8)
xlim([0,30])

imagename = sprintf('mkfigures/SINV_EPS%.2fG%d.png',EPS,Gtimes2);
saveas(gca,imagename,'png')

cd mkfigures/