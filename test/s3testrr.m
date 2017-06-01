%this code tests the calculation of 3-point correlation
clear;

%Chain structural information
N=0.1;
NRR=20;

%wavevector and structure factor
% QM=linspace(1e-5,10,500)';

% QM=linspace(1e-5,10,20)'/N;
% QMRR=linspace(1e-5,10,20)'/NRR;
QM=logspace(-3,3,20)'/N;
QMRR=logspace(-3,3,20)'/NRR;

%Chain chemical information
FA=1.00;
% ang=pi-1e-3;
ang = 2*pi/3;

%calculate s3
g3gc=zeros(length(QM),2,2,2);
g3wlc=zeros(length(QM),2,2,2);
g3rr=zeros(length(QM),2,2,2);
for ii=1:length(QM)
    ii
    Q1=QM(ii)*[1,0,0];
    Q2=transpose(rotz(ang)*Q1');
    Q3=-Q1-Q2;
    
    g3gc(ii,:,:,:)=s3gc(N,FA,Q1,Q2,Q3)/power(N,3);
    g3wlc(ii,:,:,:)=s3wlc(N,FA,Q1,Q2,Q3)/power(N,3);
    
    Q1=QMRR(ii)*[1,0,0];
    Q2=transpose(rotz(ang)*Q1');
    Q3=-Q1-Q2;
    g3rr(ii,:,:,:)=s3rr(NRR,FA,Q1,Q2,Q3)/power(NRR,3);
end

% %make plots
% figure;hold;set(gca,'fontsize',15);leg=[];
% plot(QM,g3gc(:,1,1,1).*QM*N,'k-',...
%      QM,g3gc(:,1,2,1).*QM*N,'b-',...
%      QM,g3gc(:,1,1,2).*QM*N,'r-','linewidth',2);
% plot(QM,g3wlc(:,1,1,1).*QM*N,'k--',...
%      QM,g3wlc(:,1,2,1).*QM*N,'b--',...
%      QM,g3wlc(:,1,1,2).*QM*N,'r--','linewidth',2);
% plot([0,10],[pi,pi],'k:');
% set(gca,'xscale','linear');set(gca,'yscale','linear');

% figure;hold;set(gca,'fontsize',15);
% plot(QM*N,g3wlc(:,1,1,1).*QM*N,'b-',...
%      QM*N,g3wlc(:,1,2,1).*QM*N,'b--',...
%      QM*N,g3wlc(:,1,1,2).*QM*N,'b-.','linewidth',2);
% plot(QMRR*NRR,g3rr(:,1,1,1).*QMRR*NRR,'r-',...
%      QMRR*NRR,g3rr(:,1,2,1).*QMRR*NRR,'r--',...
%      QMRR*NRR,g3rr(:,1,1,2).*QMRR*NRR,'r-.','linewidth',2);
% plot([0,50],[pi,pi],'k:');
% % axis([2,50,2,4]);
% set(gca,'xscale','linear');set(gca,'yscale','linear');
% xlabel('Lk');ylabel('SkL');box on


figure;hold;set(gca,'fontsize',15);
plot(QM*N,g3wlc(:,1,1,1),'b-')
plot(QMRR*NRR,g3rr(:,1,1,1),'r-')
plot([0,50],[pi,pi],'k:');
% axis([2,50,2,4]);
set(gca,'xscale','linear');set(gca,'yscale','linear');
xlabel('Lk');ylabel('SkL');box on