% plot Bates data with theoretical predictions
clear;
%close all

figure;
set(gca,'fontsize',20)
hold

% Collect some data
% Dataset 1: Paper from: gillard2015fluctuations
N = 39;
fa = 0.51;v = 118;b = (1-fa)*7.2+fa*8.3;
alpha = b/power(v,1/3);
Nbar = N*power(alpha,6);
CHINODT1 = 22;

% SCFT Theory
% NV = logspace(0,log10(500),100)';
NV = logspace(0,4,100)';
SCFT = ones(length(NV),1)*10.495;

% FH Theory
FH = 10.495+41.0*power(NV*alpha^6,-1/3);

% Morse Theory
MORSE = 10.495+41.0*power(NV*alpha^6,-1/3)+123*power(NV*alpha^6,-0.56);

NV = logspace(0,4,10)';
CV = power(sqrt(r2(NV)),3)./NV.*power(alpha,3);
for ii = 1:length(NV)
    [chit(ii),phase]=spinodalRG(NV(ii),CV(ii),0.5);
end
chiall = chit'.*NV;

plot(NV,SCFT,'k--')
plot(NV,FH,'b--')
plot(NV,MORSE,'r--')
plot(N,CHINODT1,'ko','markersize',10,'markerfacecolor','k');

% process figure
set(gca,'xscale','log')
% xlim([1e2,1e4])
box on