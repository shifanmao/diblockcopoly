%this code tests the calculation of 4-point correlation
clear;
addpath('../misc')
addpath(genpath('../chainstats'))

%Chain structural information
% N=100;
N=100;
NRR=20;

%wavevector and structure factor
NQ = 50;
Q0 = 1e-1;
QF = 1e2;

%QM=linspace(1e-2,2,NQ)';
QM=linspace(Q0,QF,NQ)'/N;
QMRR=linspace(Q0,QF,NQ)'/NRR;

QM = logspace(-2,2,100);
QMRR = logspace(-2,2,100)'/NRR;

%Chain chemical information
FA=1.0;
ang=pi-1e-3;
ang=pi/2;

%%%% Gaussian Chain %%%%
%calculate s4
g4gc=zeros(length(QM),2,2,2,2);
g4wlc=zeros(length(QM),2,2,2,2);
g4rr=zeros(length(QM),2,2,2,2);
for ii=1:length(QM)
    ii
    Q1=QM(ii)*[1,0,0];
    Q2=transpose(rotz(ang)*Q1');
    Q3=-Q2;
    Q4=-Q1;
    g4gc(ii,:,:,:,:)=s4gc(N,FA,Q1,Q2,Q3,Q4)./power(N,4);
%     g4wlc(ii,:,:,:,:)=s4wlc(N,FA,Q1,Q2,Q3,Q4)./power(N,4);
%     
   Q1=QMRR(ii)*[1,0,0];
   Q2=transpose(rotz(ang)*Q1');
   Q3=-Q2;
   Q4=-Q1;
    g4rr(ii,:,:,:,:)=s4rr(NRR,FA,Q1,Q2,Q3,Q4)./power(NRR,4);
end

%%make plots
%figure;hold;set(gca,'fontsize',15);
%plot(QM,g4gc(:,1,1,1,1).*QM*N,'k-',...
%     QM,g4gc(:,1,2,1,2).*QM*N,'r-','linewidth',2);
%plot(QM,g4wlc(:,1,1,1,1).*QM*N,'k--',...
%     QM,g4wlc(:,1,2,1,2).*QM*N,'r--','linewidth',2);
%axis([0,10,0,20]);
%set(gca,'xscale','linear');set(gca,'yscale','linear');
%xlabel('K');ylabel('SkL');box on
% 
% figure;hold;set(gca,'fontsize',15);
% plot(QM*N,g4gc(:,1,1,1,1).*QM*N,'k-',...
%      QM*N,g4gc(:,1,2,1,2).*QM*N,'k--','linewidth',2);
% plot(QM*N,g4wlc(:,1,1,1,1).*QM*N,'b-',...
%      QM*N,g4wlc(:,1,2,1,2).*QM*N,'b--','linewidth',2);
% %plot(QMRR*NRR,g4rr(:,1,1,1,1).*QMRR*NRR,'r-',...
% %     QMRR*NRR,g4rr(:,1,2,1,2).*QMRR*NRR,'r--','linewidth',2);
% plot([0,50],2/3*[pi,pi],'k:');
% axis([0,10,0,20]);
% set(gca,'xscale','linear');set(gca,'yscale','linear');
% xlabel('kL');ylabel('SkL');box on
% 
% 
% 
b 