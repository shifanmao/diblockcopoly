% Plot renormalized spinodal at FA=0.5 (transition to LAM phase)
% Correspond to Figure 2 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"
cd ..
addpath('functions/')

FA = 0.5;
NV = [1000,100,10,1];
CV=logspace(1,4,21);

figure;hold;set(gca,'fontsize',20)
leg=cell(length(NV),1);p=zeros(length(NV),1);
for ii = 1:length(NV)
  N = NV(ii);
  if length(NV)>1
    icol = (ii-1)/(length(NV)-1);
  else
    icol = 1;
  end
  col = [0 0 icol];

  % calculate renormalized spinodal
  [chit,phase]=spinodalRG(N,CV,FA);
  chit=reshape(chit,length(CV),1);
  chis=spinodal(N,FA);

  plot(CV.^2,ones(length(CV),1)*chis*N,'--','linewidth',2,'color',col);
  coef = (chit(1)*N-chis*N).*power(CV(1),2/3);
  
  p(ii)=plot(CV.^2,chit*N,'-','linewidth',2,'color',col);
  if N>10
      leg{ii}=strcat('\chi_S^{RG}N=',sprintf('%.2f+%.2fC^{-2/3} (N=10^%d)',chis*N,coef,log10(N)));
  else
      leg{ii}=strcat('\chi_S^{RG}N=',sprintf('%.2f+%.2fC^{-2/3} (N=%d)',chis*N,coef,N));
  end
end

set(gca,'xscale','log');box on
xlabel('C^2');ylabel('\chi_S^{RG}N');
legend(p,leg);
legend boxoff
savename = sprintf('mkfigures/figure2A.eps');
saveas(gcf,savename,'epsc')

cd mkfigures/