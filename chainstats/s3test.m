%this code tests the calculation of 3-point correlation
% clear;
addpath('misc')

%Chain structural information
NM_GS=20;

%Chain chemical information
FA=0.8;

% parameters for WLC calculations
ORDEig=12;
ORDL=12;
NumLayer=500;

%wavevector and structure factor
QM=logspace(-1,3,20)';
Q1=zeros(length(QM),1);
Q2=zeros(length(QM),1);
Q3=zeros(length(QM),1);
ang=pi;
for ii=1:length(QM)
    Q1(ii,1:3)=QM(ii)*[1,0,0];
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)');
    Q3(ii,1:3)=-Q1(ii,1:3)-Q2(ii,1:3);
end

%begin making plots
figure;hold;set(gca,'fontsize',15);leg=[];

%%%% Gaussian Chain %%%%
%calculate s3
g3=zeros(length(QM),2,2,2);
for ii=1:length(QM)
    ii
%     g3(ii,:,:,:)=s3gaussian(NM_GS,FA,Q1(ii,1:3)/NM_GS,...
%                                      Q2(ii,1:3)/NM_GS,...
%                                      Q3(ii,1:3)/NM_GS);
%     g3(ii,:,:,:)=s3wlc(NM_GS,FA,Q1(ii,1:3)/NM_GS,...
%                                 Q2(ii,1:3)/NM_GS,...
%                                 Q3(ii,1:3)/NM_GS,ORDEig,ORDL,NumLayer);
    g3(ii,:,:,:)=s3rigid(NM_GS,FA,Q1(ii,1:3)/NM_GS,...
                                Q2(ii,1:3)/NM_GS,...
                                Q3(ii,1:3)/NM_GS);
end
g3=g3/power(NM_GS,3);    

%make plots
plot(QM,g3(:,1,1,1),'k-',...
      QM,g3(:,1,2,1),'b-',...
      QM,g3(:,1,1,2),'r-','linewidth',2);
leg=[leg {['Gaussian S_{AAA}']}];
leg=[leg {['Gaussian S_{AAB}']}];

set(gca,'xscale','log');set(gca,'yscale','linear');