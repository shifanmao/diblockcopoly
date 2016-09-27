cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

FA = 0.5;
NV = logspace(0,3,20);  % number of statistical steps of total chain
PV = [1,5,10,50]*100;

figure;hold;set(gca,'fontsize',20)

% plot mean-field theory
[chis,ks,d2gam2]=spinodal(NV,FA);
plot(NV,chis.*NV,'color','k','linewidth',2);

col = 1;
for jj = 1:length(PV)
    COL = (col-1)/(length(PV)-1);
    
%     % Option 1: Keep L0 the same
%     LoverV = PV(jj);
%     CV = power(sqrt(r2(NV)),3)./NV.^3*LoverV;

%     % Option 2: Keep LP the same
%     LP3overV = PV(jj)
%     CV = power(sqrt(r2(NV)),3).*LP3overV;

    % Option 2: Keep C the same
    C = sqrt(PV(jj));
    CV = ones(length(NV),1)*C;

    chit = zeros(1,length(NV));
    for ii = 1:length(NV)
        N = NV(ii);
        C = CV(ii);
        [chit(ii),phase]=spinodalRG(N,C,0.5);
    end
    plot(NV,chit.*NV,'color',[COL 0 1-COL],'linewidth',2);
    col = col+1;
end
xlabel('N');ylabel('\chiN');

ylim([6,500]);
set(gca,'linewidth',2);box on
% legend('MF Theory','l_P^3/v = 1','l_P^3/v = 5','l_P^3/v = 10','l_P^3/v = 50')
% legend('MF Theory','L^3/Nv = 100','L^3/Nv = 500','L^3/Nv = 1000','L^3/Nv = 5000')
legend('MF Theory','C^2 = 100','C^2 = 500','C^2 = 1000','C^2 = 5000')
% title(sprintf('C^2=%.2f',C^2))

cd mkfigures/