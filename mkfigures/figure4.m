clear
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

for N = [10,50,100]
    FAV = 0.1:0.01:0.50;
    for alpha = [2, 4];
        % Renormalized phase diagrams
        [chis,chit,phase,chi13,chi36]=plotphaseRG(N,alpha,FAV);
        
        imagename = sprintf('Phase_N1e%.2f_alpha%.2f.eps',log10(N),alpha);
        filename = strcat('mkfigures/diblock-figures/figures2/',imagename);
        saveas(gca,filename,'epsc')

        % Mean-field phase diagrams
        plotphase(N,FAV);
        imagename = sprintf('MFPhase_N1e%.2f.eps',log10(N));
        filename = strcat('mkfigures/diblock-figures/figures2/',imagename);
        saveas(gca,filename,'epsc')
    end
end

cd mkfigures/