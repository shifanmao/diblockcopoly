clear;

% Figure 1
filename = 'fig1rrplot';
outfile = fopen(filename, 'wt');

NV = logspace(-1,3,81);

% alphaV = [1,2,4,8,16];
alphaV = 1;
for alpha = alphaV
    for nn = 1:length(NV)
        N = NV(nn);
        [chit,phase]=spinodalrrRG(N,alpha,0.5);
        
        result = [alpha, N, chit*N];
        fprintf(outfile,'%.2f, %.4f, %.4f\n',result);
    end
end

% filename = 'fig1plot';
% fig1 = load(filename);
% 
% figure;hold
% alphaV = 1:5;
% for alpha = alphaV
%     NV = fig1(fig1(:,1)==alpha, 2);
%     CHITN = fig1(fig1(:,1)==alpha, 3);
%     plot(NV, CHITN)
% end
% set(gca,'xscale','log')
% axis([1,1000,5,30])

% % Figure 2
% filename = 'fig3plot';
% outfile = fopen(filename, 'wt');
% 
% NV = [10, 50, 100];
% FAV = linspace(.1, .9, 41);