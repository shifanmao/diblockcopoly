% clear;close all
% 
% result = [];
% NV = logspace(-1,4,21);
% for N = NV
%     % Figure 5: renormalized spinodal at FA=0.5
%     CV=logspace(1,4,21)';
%     [chit,phase]=spinodalRG(N,CV,0.5);
%     chit=reshape(chit,length(CV),1);
%     
%     result = [result; repmat(log10(N),length(CV),1),log10(CV),chit*N];
% end
% dlmwrite('data/odt',result,'delimiter','\t','precision',4);


% make a plot
data = load('data/odt');
NV = logspace(-1,4,21);
CV = logspace(1,4,21)';
CHITN = zeros(21,21);
for ii = 1:21
	N = NV(ii);
    for jj = 1:21
        C = CV(jj);
%         ind = find( abs(data(:,1)-log10(N))<1e-2 & abs(data(:,2)==log10(C))<1e-2);
        ind = find( abs(data(:,1)-log10(N))<1e-2 & data(:,2)==log10(C));
        CHITN(jj,ii) = data(ind,3);
    end
end

figure;set(gca,'fontsize',50)
surf(NV,CV.^2,CHITN,'edgecolor','none','LineStyle','none','FaceLighting','phong');
set(gca,'xscale','linear');set(gca,'yscale','linear');
%    xlim([min(TV),max(TV)]);ylim([min(QM),max(QM)]);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('N');
ylabel('C^2');
zlabel('\chi_tN')
% colorbar;
view([0,90]);colormap(jet)
% set(gca, 'CLim', [0,2.]);