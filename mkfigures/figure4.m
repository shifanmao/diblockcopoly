clear
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

for N = [10,50,100]
    FAV = linspace(.1,.5,41);
    for twoLPoverV = [2,4];
        C = power(sqrt(r2(N)),3)/N*power(twoLPoverV,3);
        [chis,chit,phase,chi13,chi36]=plotphaseRG(N,C,FAV);
        ylim([8,30])
        imagename = sprintf('Phase_N1e%.2f_LPoV%.2f.eps',log10(N),twoLPoverV);
        filename = strcat('mkfigures/figures2/',imagename);
        saveas(gca,filename,'epsc')

        plotphase(N,FAV);
        ylim([8,30])
        imagename = sprintf('MFPhase_N1e%.2f.eps',log10(N));
        filename = strcat('mkfigures/figures2/',imagename);
        saveas(gca,filename,'epsc')
    end
end

cd mkfigures/