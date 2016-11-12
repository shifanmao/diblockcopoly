clear
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

% for N = [100,562.34,1e3]
for N = [562.34,1e3]
FAV = linspace(.1,.5,41);
for twoLPoverV = [4];

C = power(sqrt(r2(N)),3)/N*power(twoLPoverV,3);
[chis,chit,phase,chi13,chi36]=plotphaseRG(N,C,FAV);
imagename = sprintf('Phase_N1e%.2f_LPoV%.2f.eps',log10(N),twoLPoverV);
filename = strcat('mkfigures/figures2/',imagename);
saveas(gca,filename,'epsc')

% filename = '';
% plotphase(N,FAV);
% imagename = sprintf('MFPhase_N1e%.2f.eps',log10(N));
% filename = strcat('mkfigures/figures2/',imagename);
% saveas(gca,filename,'epsc')
end
end

cd mkfigures/