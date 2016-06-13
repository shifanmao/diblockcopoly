function R2=r2(NM)
% Mean-squared end-to-end distance of wormlike chain
% Usage :: [r2]=r2wlc(N)
% Output :: r2 = mean-squared end-to-end distance of wormlike chains
%           in the limit NM>=1e4, use Gussian chain model
%           in the limit NM<=1e-4, use perfectly rigid rod
% Input :: N = number of Kuhn steps of wormlike chain

if NM>=1e4  % Gaussian chain limit
    R2 = NM;
elseif NM<=1e-4  % Rigid rod limit
    R2 = NM^2;
else  % Worm-like chain
    R2 = NM-(1/2)*(1-exp(-2*NM));
end
