clear
close all

d=3.0;
K0=1e-3;
KF=1e5;
NK=1000;
K=transpose(logspace(log10(K0),log10(KF),NK));

L1 = 2;
L2 = 1;
M = 1;

ORDL=3;
ORDM=3;
ORD=20;
ResLayer=500;

R=zeros(NK,ORD);
G=zeros(NK,ORD);
GA=zeros(NK,ORD);

for I=1:NK
%    EI=MatRoots(K(I),d,ORD);
%    R(I,:)=transpose(EI);

    EI=Eigenvalues(K(I),ORD,ORDM);
    R(I,:)=transpose(EI(:,M+1));
end

for I=1:NK
    RESI=Residues(K(I),R(I,:),ORD,ORDL,ResLayer,ORDM);
    G(I,:)=RESI(L1+1,L2+1,M+1,:);
end


for I=1:ORD
    COL=(I-1)/(ORD-1);
    figure(1)
    semilogx(K,real(R(:,I)),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    figure(2)
    semilogx(K,imag(R(:,I)),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    
    figure(3)
    loglog(K,abs(G(:,I)),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
        
    figure(4)
    loglog(K,abs(imag(G(:,I))),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on    
end

%figure(3)
%axis([1e-3 1e5 0 1])
