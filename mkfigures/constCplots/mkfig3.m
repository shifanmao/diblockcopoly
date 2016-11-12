% Make the phase diagrams at different chain rigidity and chain effective densities
% Correspond to Figure 3 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"

cd ..
addpath('functions/')

NV = [1,10,100,1000];
CV = [100,1000];
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition

for N = NV
  plotphase(N,FAV);
  set(gca,'Xtick',0.1:0.1:0.9)
  set(gca,'XtickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})

  savename = sprintf('mkfigures/N%dCinf.eps',N);
  saveas(gcf,savename,'epsc')
end

for N = NV
  for C = CV
    plotphaseRG(N,C,FAV);
    set(gca,'Xtick',0.1:0.1:0.9)
    set(gca,'XtickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})

    savename = sprintf('mkfigures/N%dC%d.eps',N,C);
    saveas(gcf,savename,'epsc')
  end
end

cd mkfigures