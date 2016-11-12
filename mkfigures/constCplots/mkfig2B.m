% Plot renormalized spinodal at FA=0.5 (transition to LAM phase)
% Correspond to Figure 2 in the manuscript "Diblock Phase Behavior: chain semiflexibility and density fluctuation effects"
cd ..
addpath('functions/')
addpath(genpath('chainstats/'))
addpath('misc/')

% Figure 2B: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(-1,3,51);  % number of statistical steps of total chain
[chis,ks,d2gam2]=spinodal(NV,0.5);

figure;set(gca,'fontsize',20);hold
%plot(NV,chis.*NV,'k-','linewidth',2);
% vary color
x = NV;
y = chis.*NV;
z = zeros(size(x));
col = x;  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2,'linewidth',3);
map = zeros(100,3);
map(:,3) = logspace(0,-9,100)';
colormap(map)
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2,'linewidth',3);

set(gca,'xscale','log');box on;
xlim([1e-1,1e3]);ylim([6,11.5])
xlabel('N');ylabel('\chi_S^{MF}N')
plot(NV,ones(length(NV),1)*6.135,'b--','linewidth',2)
plot(NV,ones(length(NV),1)*10.495,'k--','linewidth',2)
% text(40,6.135+0.31,'(\chi_S^{MF}N)_{RR}=6.135','FontSize',20,'color','b')
% text(0.15,10.495+0.31,'(\chi_S^{MF}N)_{GC}=10.495','FontSize',20,'color','k')
text(2.5,6.135-0.91,'(\chi_S^{MF}N)_{RR}=6.135',...
    'interpreter','tex','FontSize',20,'color','b')
text(2.5,10.495+0.31,'(\chi_S^{MF}N)_{GC}=10.495','FontSize',20,'color','k')


savename = sprintf('mkfigures/figure2B.eps');
saveas(gcf,savename,'epsc')

% Figure 2C: renormalized spinodal coefficient at FA=0.5
FA = 0.5;
NV = logspace(3,0,81);
CV = 10;

leg=cell(length(NV),1);p=zeros(length(NV),1);
coef=zeros(length(NV),1);
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

  coef(ii) = (chit(1)*N-chis*N).*power(CV(1),2/3);
end

figure;set(gca,'fontsize',20);hold
semilogx(NV,coef,'k-','linewidth',2);
% plot empirical solutions
plot(NV,ones(length(NV),1)*41.01,'k--','linewidth',2)
plot(logspace(0,0.5,3),1.05*74.16*power(logspace(0,0.5,3),-0.27),'k--','linewidth',2)
text(2,70,'f(N)~N^{-0.27}','FontSize',20)
text(150,43,'f(N)=41.01','FontSize',20)
set(gca,'xscale','log');set(gca,'yscale','log');box on
xlabel('N');ylabel('f(N)');
ylim([40,80]);
savename = sprintf('mkfigures/figure2C.eps');
saveas(gcf,savename,'epsc')

cd mkfigures/