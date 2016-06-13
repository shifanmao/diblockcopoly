function s2=s2gc(NM,FA,k)
% Calculate the Fourier transform of the Green function
% for the gaussian chain in d-dimension
%
% Andrew Spakowitz (4/14/15)

d=3;    %dimension
FB=1-FA;
R =-NM*k*k/(2*d);

s2=zeros(2,2);
if abs(R)<1e-10
    s2(1,1)=FA*FA;
    s2(2,2)=FB*FB;
    s2(1,2)=FA*FB;
else
% Case 1 :: A1==A2
    % on same monomer (S integrals)
    valeqA = R.^(-2).*(-1+exp(FA.*R)-FA.*R);
    valeqB = R.^(-2).*(-1+exp(FB.*R)-FB.*R);

    s2(1,1)=s2(1,1)+2*valeqA*(NM^2);
    s2(2,2)=s2(2,2)+2*valeqB*(NM^2);

% Case 2 :: A1~=A2
    % on same monomer (S integrals)
    valeqAB = (1+exp(R)+(-1).*exp(FA.*R)-exp(R-FA.*R)).*R.^(-2);
    % J1<J2
    s2(1,2)=s2(1,2)+valeqAB*(NM^2);
end
s2(2,1)=s2(1,2);

end