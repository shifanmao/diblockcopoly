%this code tests the calculation of 4-point correlation
clear;
addpath('../misc')
addpath('integrals/')
addpath('eigcalc/')

%Chain structural information
NM=20;

%Chain chemical information
FA=1.0;

% parameters for WLC calculations
ORDEig=2;
ORDL=1;
NumLayer=500;

%wavevector and structure factor
QM=logspace(-1,2,20)';
Q1=zeros(length(QM),1);
Q2=zeros(length(QM),1);
Q3=zeros(length(QM),1);
Q4=zeros(length(QM),1);
ang=pi;
for ii=1:length(QM)
    Q1(ii,1:3)=QM(ii)*[1,0,0];
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)');
    Q3(ii,1:3)=-Q2(ii,1:3);
    Q4(ii,1:3)=-Q1(ii,1:3);
end

%begin making plots
figure;hold;set(gca,'fontsize',15);leg=[];

%%%% Gaussian Chain %%%%
%calculate s4
g4=zeros(length(QM),2,2,2,2);
for ii=1:length(QM)
   ii
%    g4(ii,:,:,:,:)=s4gc(NM,FA,Q1(ii,1:3),Q2(ii,1:3),Q3(ii,1:3),Q4(ii,1:3));
%    g4(ii,:,:,:,:)=s4wlc(NM,FA,Q1(ii,1:3),Q2(ii,1:3),Q3(ii,1:3),Q4(ii,1:3));
                                  %,ORDEig,ORDL,NumLayer);
   g4(ii,:,:,:,:)=s4rr(NM,FA,Q1(ii,1:3),Q2(ii,1:3),Q3(ii,1:3),Q4(ii,1:3));
end
g4=g4./power(NM,4);

%make plots
plot(QM,g4(:,1,1,1,1),'kx',...
     QM,g4(:,1,2,1,2),'r-','linewidth',2);
leg=[leg {'Gaussian'}];

legend(leg);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.3]);
xlabel('Normalized Wavevector kL');
ylabel('Normalized Structure Factor S/L^4')