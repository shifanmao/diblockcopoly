clear
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

for N = [10,100]
    FAV = linspace(.1,.5,41);
    for alpha = [4];
        [chis,chit,phase,chi13,chi36]=plotphaseRG(N,alpha,FAV);
%         ylim([8,30])
        ylim([5,20])
        imagename = sprintf('Phase_N1e%.2f_alpha%.2f.eps',log10(N),alpha);
        filename = strcat('mkfigures/figures2/',imagename);
        saveas(gca,filename,'epsc')

        plotphase(N,FAV);
%         ylim([8,30])
        ylim([5,20])
        imagename = sprintf('MFPhase_N1e%.2f.eps',log10(N));
        filename = strcat('mkfigures/figures2/',imagename);
        saveas(gca,filename,'epsc')
    end
end

cd mkfigures/