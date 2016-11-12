clear all
cd ../

addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')


%%%%%%%%%% MAKE A PLOT %%%%%%%%%%
NbarV = logspace(2,8,100);
chiMFV = 10.495 + zeros(length(NbarV),1);
chiFHV = 10.495 + 41.0*power(NbarV,-1/3);
chiGLV = 10.495 + 41.0*power(NbarV,-1/3) + 123.0*power(NbarV,-0.56);

figure;hold;set(gca,'fontsize',20)
plot(NbarV,chiMFV,'k--')
plot(NbarV,chiFHV,'k-.')
plot(NbarV,chiGLV,'k-')
set(gca,'xscale','log')
xlabel('Nbar');
ylabel('\chiN_{ODT}')

% load Bates Data
NbarData = [2.77e4,4270,3250];
NData    = [1925,791,427];
EData    = [0.2,2.0,0.3];
errorbar(NbarData,chiODT,EData,'ro')

% start calculation
CHIODT   = zeros(1,length(NbarData));
for ii = 1:length(NbarData)
    Nbar = NbarData(ii);
    N = NData(ii);
    
    % ODTs reported by Bates
    chiODT(ii) = 10.90+41.0*power(Nbar,-1/3);

    %%%%%%%%%% PREDICT PHASE TRANSITIONS %%%%%%%%%%

    % F-H Prediction
    chiFH = 10.495 + 41.0*power(Nbar,-1/3);
    % Glaser Predictions
    chiGL = 10.495 + 41.0*power(Nbar,-1/3) + 123.0*power(Nbar,-0.56);

%     plot(Nbar,chiODT(ii),'ro','markerfacecolor','r');
%     plot(Nbar,chiFH,'ko')
%     plot(Nbar,chiGL,'ko')
    
    % My Prediction
    LPoverV = power(Nbar/N,1/6)
    C = power(sqrt(r2(N)),3)./N.*power(LPoverV,3);
    [chit,phase]=spinodalRG(N,C,0.5);
    chiMY = chit*N;
    plot(Nbar,chiMY,'bo','markerfacecolor','b');
end


















% %%%%%%%%%% MAKE A PLOT %%%%%%%%%%
% NbarV = logspace(2,8,100);
% chiMFV = 10.495 + zeros(length(NbarV),1);
% chiFHV = 10.495 + 41.0*power(NbarV,-1/3);
% chiGLV = 10.495 + 41.0*power(NbarV,-1/3) + 123.0*power(NbarV,-0.56);
% 
% figure;hold;set(gca,'fontsize',20)
% plot(NbarV,chiMFV,'k--')
% plot(NbarV,chiFHV,'k-.')
% plot(NbarV,chiGLV,'k-')
% set(gca,'xscale','log')
% xlabel('Nbar');
% ylabel('\chiN_{ODT}')
% 
% % load Bates Data
% NbarData = [2.77e4,4270,3250];
% NData    = [1925,791,427];
% EData    = [0.2,2.0,0.3];
% CHIODT   = zeros(1,length(NbarData));
% 
% for ii = 1:length(NbarData)
%     Nbar = NbarData(ii);
%     N = NData(ii);
%     
%     % ODTs reported by Bates
%     chiODT(ii) = 10.90+41.0*power(Nbar,-1/3);
% 
%     %%%%%%%%%% PREDICT PHASE TRANSITIONS %%%%%%%%%%
% 
%     % F-H Prediction
%     chiFH = 10.495 + 41.0*power(Nbar,-1/3);
%     % Glaser Predictions
%     chiGL = 10.495 + 41.0*power(Nbar,-1/3) + 123.0*power(Nbar,-0.56);
% 
% %     plot(Nbar,chiODT(ii),'ro','markerfacecolor','r');
% %     plot(Nbar,chiFH,'ko')
% %     plot(Nbar,chiGL,'ko')
%     
%     % My Prediction
%     LPoverV = power(Nbar/N,1/6)
%     C = power(sqrt(r2(N)),3)./N.*power(LPoverV,3);
%     [chit,phase]=spinodalRG(N,C,0.5);
%     chiMY = chit*N;
%     plot(Nbar,chiMY,'bo','markerfacecolor','b');
% end
% errorbar(NbarData,chiODT,EData,'ro')

cd test/