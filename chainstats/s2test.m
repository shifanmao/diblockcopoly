%this code tests the calculation of 2-point correlation
clear;

%Check ABs (if chkab==1, check combinations of SAB, SBA, etc.)
chkab=10;

%Chain structural information
NM_GS=50;

%Chain chemical information
FA=0.8;

% parameters for WLC calculations
ORDEig=10;
NumLayer=500;

%wavevector
k=logspace(-1,3,50)';

%begin making plots
figure;hold;set(gca,'fontsize',15);leg=[];

%%%% Gaussian Chain %%%%
%calculate s2
g2=zeros(length(k),2,2);
for ii=1:length(k)
    ii
    g2(ii,:,:) = s2gc(NM_GS,FA,k(ii)/NM_GS);
%     g2(ii,:,:) = s2wlc(NM_GS,FA,k(ii)/NM_GS,ORDEig,NumLayer);
%     g2(ii,:,:) = s2rigid(NM_GS,FA,k(ii)/NM_GS);
end
g2 = g2./power(NM_GS,2);

leg=[];
%make plots
plot(k,g2(:,1,1),'k-',k,g2(:,1,2),'k--','linewidth',2);
leg=[leg {'Gaussian S_{AA}'}];
leg=[leg {'Gaussian S_{AB}'}];

legend(leg);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.2]);
xlabel('Normalized Wavevector kL');ylabel('Normalized Structure Factor S/L^2')