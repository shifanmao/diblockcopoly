cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

FA = 0.5;
NV=logspace(0,3,21);
PV = [2,5,20,50]*1e5;

figure;hold;set(gca,'fontsize',20)

% plot mean-field theory
[chis,ks,d2gam2]=spinodal(NV,FA);
plot(NV,chis.*NV,'color','k','linewidth',2);

col = 1;
leg = {};
for jj = 1:length(PV)
    COL = (col-1)/(length(PV)-1);

    % Keep LP (aspect ratio) the same
    LP3overV = PV(jj)
    CV = power(sqrt(r2(NV)),3)./N.*LP3overV;

    chit = zeros(1,length(NV));
    for ii = 1:length(NV)
        N = NV(ii);
        C = CV(ii);
        [chit(ii),phase]=spinodalRG(N,C,0.5);
    end
    plot(NV,chit.*NV,'color',[COL 0 1-COL],'linewidth',2);
    col = col+1;
    
    leg{jj} = strcat('(2l_P)^3/v=',sprintf('%d',PV(jj)));
end
xlabel('N');ylabel('\chiN');

% % plot simulation results
% C2sim = 1525.9;
% SIMEPS = 32*[0.05,0.1,0.5,1.0];
% SIMODT1 = [16.361,15.958,15.914,16.094];
% SIMODT2 = [16.672,16.533,16.311,16.374];
% SIMODT = (SIMODT1+SIMODT2)/2;
% COL = (3-1)/(length(PV)-1);
% plot(SIMEPS,SIMODT,'--','linewidth',2,'color',[COL 0 1-COL]);
% plot(SIMEPS,SIMODT1,'s','linewidth',2,'color',[COL 0 1-COL]);
% plot(SIMEPS,SIMODT2,'s','linewidth',2,'color',[COL 0 1-COL]);

ylim([6,500]);
set(gca,'linewidth',2);box on
legend(leg);
cd mkfigures/