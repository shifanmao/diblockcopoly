function s2inv=s2inverse(N,FA,k)
%% calculates the inverse of pair correlation function in Fourier space
% usage s2inv=s2inverse(N,FA,k)
% Return:
%    val, inverse s2
% Parameters:
%    N, number of Kuhn steps per monomer
%    FA, fraction of A type monomers
%    k, magnitude of wavevector in unit of 1/contour length
%       k input can be a vector

s2inv=zeros(2,2);
MIN=1e-4;
NRR=50;  % rigid rod bead discretization

% Calculate the s matrix
if N>=1e4  % Gaussian chain limit
    s2=s2gc(N,FA,k);
elseif N<=1e-2  % Rigid rod limit
    s2=s2rr(NRR,FA,k);
else
    s2=s2wlc(N,FA,k);
end
DET=s2(1,1,:)*s2(2,2,:)-s2(1,2,:)*s2(2,1,:);

if abs(k)<MIN
    if N<=1e-2
      s2inv=ones(2,2)./power(NRR,2);
    else
      s2inv=ones(2,2)./power(N,2);
    end
else
    s2inv(1,1) = s2(2,2)./DET;
    s2inv(1,2) = -s2(1,2)./DET;
    s2inv(2,1) = -s2(2,1)./DET;
    s2inv(2,2) = s2(1,1)./DET;
end
