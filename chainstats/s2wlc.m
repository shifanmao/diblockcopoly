function s2=s2wlc(NM,FA,k)
%% function s2wlc :: Calculate the Fourier transform of the Green function
% for the worm-like chain in d-dimension
% Usage: s2=s2wlc(NM,FA,k)
% Andrew Spakowitz (4/14/15)

% parameters for worm-like chain calculations
% d=3;  % number of dimensions
ORDmax=20;  % maximum number of eigenvalues
ResLayer=500;  % number of residual layers
FB=1-FA;

% calculate the eigenvalues
R=Eigenvalues(k,ORDmax,1);
NR=ORDmax;

% get the residues for all roots of each k(j)
Residue=Residues(k,R(1:NR),NR,1,ResLayer,1);

s2=zeros(2,2);
if abs(R(1))<1e-10
    s2(1,1)=FA*FA;
    s2(2,2)=FB*FB;
    s2(1,2)=FA*FB;
else
    [j0,dj0]=CalRes0(k,ResLayer);
    G0 = 1/j0;
    pG0 = -dj0/j0^2;
    
    s2(1,1)=s2(1,1) + 2*(NM*FA*G0 + pG0);
    s2(2,2)=s2(2,2) + 2*(NM*FB*G0 + pG0);
    s2(1,2)=s2(1,2) - pG0;
    
    for I=1:NR
        R(I)=R(I)*NM;

        % Case 1 :: A1==A2
        valAA = R(I).^(-2).*exp(FA.*R(I));
        valBB = R(I).^(-2).*exp(FB.*R(I));

        s2(1,1)=s2(1,1)+2*Residue(I)*valAA*(NM^2);
        s2(2,2)=s2(2,2)+2*Residue(I)*valBB*(NM^2);

        % Case 2 :: A1~=A2
        valAB = (exp(R(I))-exp(FA.*R(I))-exp(R(I)-FA.*R(I))).*R(I).^(-2);
        s2(1,2)=s2(1,2)+Residue(I)*valAB*(NM^2);
    end
end
s2(2,1)=s2(1,2);

s2(imag(s2) < 1e-2) = real(s2);

end

function [j0,dj0]=CalRes0(k,ResLayer)
d = 3;

if k<=1
    
    % residual using continued fraction for small k
    p=0;
    W=p+(ResLayer+d-2)*ResLayer;
    Wprime=1;
    for L=ResLayer:-1:1
        AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
        Wprime=1-AL^2*Wprime/W^2;
        PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
        W=PLm+AL^2/W;
    end
    dj0=Wprime;
    j0=W;
else
    
    % residual using continued fraction for large k
    p=0;
    W=(p+(ResLayer+d-2)*ResLayer)/k;
    Wprime=1/k;
    for L=ResLayer:-1:1
        AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
        Wprime=1/k-AL^2*Wprime/W^2;
        PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
        W=PLm/k+AL^2/W;
    end
    dj0=k*Wprime;
    j0=k*W;
end
end