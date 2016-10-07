N1 = 1;
N2 = 1000;
CHI = linspace(-1,15,100);

figure;set(gca,'fontsize',20);hold
[chis,ks,d2gam2]=spinodal(N1,0.5);
plot(CHI,-2*CHI+2*chis*N1,'b','linewidth',2)
plot(chis*N1,0,'bo','markerfacecolor','b','markersize',8)

[chis,ks,d2gam2]=spinodal(N2,0.5);
plot(CHI,-2*CHI+2*chis*N2,'k','linewidth',2)
plot(chis*N2,0,'ko','markerfacecolor','k','markersize',10)

xlim([0,12])
ylim([0,25])
xlabel('\chiN')
ylabel('N<\psi^2>_{MF}^{-1}')
box on