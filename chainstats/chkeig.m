clear
close all
addpath('eigcalc')

d=3.0;
K0=1e-2;
KF=1e2;
NK=1000;
K=transpose(logspace(log10(K0),log10(KF),NK));

ORD=10;
ResLayer=500;

ORDL = 1;
ORDMU = 1;

R=zeros(NK,ORD);
smallK=zeros(NK,ORD);
largeK=zeros(NK,ORD);
largeKV=zeros(NK,ORD,1,ORD);
G=zeros(NK,ORD);

for I=1:NK    
    EI=Eigenvalues(K(I),ORD,ORDL);
    R(I,:)=transpose(EI(:,ORDL));
end

for I=1:NK
    EigK = transpose(R(I,:));
    [smallKV,largeKV]=Residues_two(K(I),EigK,ORD,ORDL,ResLayer,ORDMU);
    smallK(I,:)=smallKV(1,1,1,:);
    largeK(I,:)=largeKV(1,1,1,:);
end

for I=ORDL:ORD
    COL=(I-1)/(ORD-1);
    figure(1)
    semilogx(K,real(R(:,I)),'.-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    figure(2)
    semilogx(K,imag(R(:,I)),'.-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on

    % find cross-over K
    x = largeK(:,I);
    d2x = diff(diff(log(x)));
    kc1 = K(find(abs(d2x)>10,1));  % remember to include k diff.
    if ~isempty(kc1)
      indc1 = find(K>kc1*10^(0.1),1);
    else
      indc1 = 1;
    end

    y = abs(largeK(:,I)-smallK(:,I));
    kc2 = K(find(y==min(y)));
    if ~isempty(kc2)
      indc2 = find(K>kc2*10^(0.1),1);
    else
      indc2 = 1;
    end
    ind = max([indc1,indc2]);
    kc = K(ind);
    G(:,I) = [smallK(1:ind,I);largeK(ind+1:end,I)];

    figure(3)
    loglog(K,abs(real(G(:,I))),'.-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    figure(4)
    loglog(K,abs(imag(G(:,I))),'.-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on

%    figure(3)
%    loglog(K,abs(real(largeK(:,I))),'-','LineWidth',2,'Color',[COL 0 1-COL])
%    hold on
%    figure(4)
%    loglog(K,abs(imag(largeK(:,I))),'-','LineWidth',2,'Color',[COL 0 1-COL])
%    hold on

%    figure(3)
%    loglog(K,abs(real(smallK(:,I))),'x','LineWidth',2,'Color',[COL 0 1-COL])
%    hold on
%    figure(4)
%    loglog(K,abs(imag(smallK(:,I))),'x','LineWidth',2,'Color',[COL 0 1-COL])
%    hold on
end
