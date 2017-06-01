function plotphase_wsolvent(FAVV, PHIPV, EIGV, EIG, KSV, N)

arrow_length = 0.025;

hold on
plotvectors(FAVV, PHIPV, EIGV(:, :, :), arrow_length);
plotvectors(FAVV, PHIPV, -EIGV(:, :, :), arrow_length);

plotsurf(FAVV, PHIPV, 1./EIG/N);
axis([-.1,2/sqrt(3)+.1,-.1,1.1])

box on
axis off
grid off
set(gca,'fontsize',18)

colormap jet
h = colorbar;
% ylabel('h, $\langle \bf{\psi}^{\bf{\top}}(\vec{q}) \bf{\psi}(-\vec{q}) \rangle/Nv$',...
%             Interpreter', 'Latex')
caxis([0,0.5])
% dif=0.05;axis([0-dif,2/sqrt(3)+dif, -dif, 1+dif])
% plot3(1/sqrt(3), 0.5, 1e6, 'ko', 'markerfacecolor', 'w', 'markersize', 10)
% plot3(1.25/sqrt(3), 0.5, 1e6, 'ko', 'markerfacecolor', 'w', 'markersize', 10)

hold off