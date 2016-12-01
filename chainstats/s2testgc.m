%this code tests the calculation of 2-point correlation
clear;

%Chain structural information
N=100;
NRR=20;

%Chain chemical information
FA=1.0;

% parameters for WLC calculations
ORDEig=10;
NumLayer=500;

%wavevector
k=linspace(1e-5,20,500)';

%%%% Gaussian Chain %%%%
%calculate s2
g2gc=zeros(length(k),2,2);
g2wlc=zeros(length(k),2,2);
g2rr=zeros(length(k),2,2);
for ii=1:length(k)
    g2gc(ii,:,:) = s2gc(N,FA,k(ii))/power(N,2);
    g2wlc(ii,:,:) = s2wlc(N,FA,k(ii))/power(N,2);
%     g2rr(ii,:,:) = s2rr(NRR,FA,k(ii))/power(NRR,2);
end

%make plots
figure;hold;set(gca,'fontsize',15);
plot(k,g2gc(:,1,1).*k*N,'k-',...
     k,g2gc(:,1,2).*k*N,'k--','linewidth',2);
plot(k,g2wlc(:,1,1).*k*N,'b-',...
     k,g2wlc(:,1,2).*k*N,'b--','linewidth',2);
plot([0,10],[pi,pi],'k:');
set(gca,'xscale','linear');set(gca,'yscale','linear');
axis([0,10,0,20]);
xlabel('K');ylabel('SkL')
% 
% figure;hold;set(gca,'fontsize',15);
% plot(k*N,g2wlc(:,1,1).*k*N,'b-',...
%      k*N,g2wlc(:,1,2).*k*N,'b--','linewidth',2);
% plot(k*NRR,g2rr(:,1,1).*k*NRR,'r-',...
%      k*NRR,g2rr(:,1,2).*k*NRR,'r--','linewidth',2);
% plot([0,20],[pi,pi],'k-.');
% axis([2,20,2,7]);
% xlabel('K');ylabel('SkL')