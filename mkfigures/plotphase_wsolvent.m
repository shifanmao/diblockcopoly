function plotphase_wsolvent(N, KV, CHI, IEIG)

FAV = linspace(.1, .9, 21);
PHIPV = linspace(.001, .999, 51);
EIGV = zeros(length(FAV), length(PHIPV), 2);
EIG = zeros(length(FAV), length(PHIPV));

% filename='data/newgamdata';
% outfile = fopen(filename, 'wt');
for jj = 1:length(PHIPV)
    PHIP = PHIPV(jj)
    for ii = 1:length(FAV)
        FA = FAV(ii);

        [KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, KV, PHIP, CHI, IEIG);
        EIGV(ii, jj, 1:2) = EIGVS;
        EIG(ii, jj) = EIGS;
    end
end


filename = sprintf('../data/N%.2fCHIABN%.2fIEIG%d.mat', N, CHI(1,2)*N, IEIG);
save(filename, 'FAV', 'PHIPV', 'EIGV', 'EIG');

plotvectors(FAV, 1-PHIPV, EIGV*0.01)

hold on

plotvectors(FAV, 1-PHIPV, -EIGV*0.01)
plotsurf(FAV, 1-PHIPV, 1./EIG/N)
axis([0,2/sqrt(3),0,1])
if IEIG == 1
    caxis([0,.2])
else
    caxis([0,1e-4])
%     caxis([0,1e-2])
end

if CHI(1,1) > 0
    title(strcat('\chiN', sprintf('=%.2f', CHI(1,1)*N)))
else
    title(strcat('\chiN', sprintf('=%.2f', CHI(1,2)*N)))
end

box on
colormap jet
colorbar
axis off
grid off
set(gca,'fontsize',18)

hold off