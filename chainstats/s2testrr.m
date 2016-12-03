%this code tests the calculation of 2-point correlation
clear;

%Chain structural information
N=1;
NRR=20;

%wavevector
k=linspace(1e-5,20,100)'/N;
kRR=linspace(1e-5,20,100)'/NRR;

%Chain chemical information
FA=1.0;

%%%% Gaussian Chain %%%%
%calculate s2
g2gc=zeros(length(k),2,2);
g2wlc=zeros(length(k),2,2);
g2rr=zeros(length(k),2,2);
for ii=1:length(k)
%     g2gc(ii,:,:) = s2gc(N,FA,k(ii))/power(N,2);
    g2wlc(ii,:,:) = s2wlc(N,FA,k(ii))/power(N,2);
    g2rr(ii,:,:) = s2rr(NRR,FA,kRR(ii))/power(NRR,2);
end

%make plots
% figure;hold;set(gca,'fontsize',15);
% plot(k,g2gc(:,1,1).*k*N,'k-',...
%      k,g2gc(:,1,2).*k*N,'k--','linewidth',2);
% plot(k,g2wlc(:,1,1).*k*N,'b-',...
%      k,g2wlc(:,1,2).*k*N,'b--','linewidth',2);
% plot([0,10],[pi,pi],'k:');
% set(gca,'xscale','linear');set(gca,'yscale','linear');
% axis([0,10,0,20]);
% xlabel('K');ylabel('SkL')

figure;hold;set(gca,'fontsize',15);
plot(k*N,g2wlc(:,1,1).*k*N,'b-',...
     k*N,g2wlc(:,1,2).*k*N,'b--','linewidth',2);
plot(kRR*NRR,g2rr(:,1,1).*kRR*NRR,'r-',...
     kRR*NRR,g2rr(:,1,2).*kRR*NRR,'r--','linewidth',2);
plot([0,20],[pi,pi],'k:');
axis([2,20,2,3.5]);
xlabel('Lk');ylabel('SkL');box on