clear
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')
figure;hold;set(gca,'fontsize',20)

FA = 0.5;
NV=logspace(-1,6,101);

% Keep LP (aspect ratio) the same
PV = fliplr([1,2,3,4,5].^3);

% plot mean-field theory
[chis,ks,d2gam2]=spinodal(NV,FA);
plot(NV,chis.*NV,'color','k','linewidth',2);
plot(NV,ones(length(NV),1)*10.45,'k--','linewidth',2);

col = 1;
leg = {};
for jj = 1:length(PV)
    COL = (col-1)/(length(PV)-1);
    chit = zeros(1,length(NV));

    % Keep LP (aspect ratio) the same
    LP3overV = PV(jj)
    
    for ii = 1:length(NV)
        N = NV(ii);
        C = power(sqrt(r2(N)),3)./N.*LP3overV;
        [chit(ii),phase]=spinodalRG(N,C,0.5);
    end
    plot(NV,chit.*NV,'color',[COL 0 1-COL],'linewidth',2);
    col = col+1;
    leg{jj} = strcat('(2l_P)^3/V=',sprintf('%d',PV(jj)));
    
    % calculate Fredrickson-Helfand Theory
    chiNFH = 10.49+41.0*power(LP3overV^2*NV,-1/3);
    plot(NV,chiNFH,'--','color',[COL 0 1-COL],'linewidth',2)

    result = [NV',chit'.*NV'];
    filename = strcat('data/chiODT',sprintf('_LP3overV%.2f',LP3overV));
    dlmwrite(filename,result,'precision','%.3f');
end
xlabel('N');ylabel('\chi_{ODT}^{OL}N');

ylim([5,25])
set(gca,'xscale','log')
set(gca,'linewidth',2);box on
% legend(leg);
cd mkfigures/