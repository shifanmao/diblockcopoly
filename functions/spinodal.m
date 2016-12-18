function [chis,ks,d2gam2]=spinodal(N,FAV)
%% spinodal.m :: This function predicts diblock copolymer phase transition spinodal
% and critical wavelength of quadratic density fluctuation in the well-mixed state
% Usage: [chis,ks,d2gam2]=spinodal(N,FAV)
% Inputs:
%   N, number of statistical steps of total chain
%   FAV, range of A-type monomer fractions
% Outputs:
%   chis, Flory-Huggins parameter at spinodal
%   ks, critical wavelength of quadratic density fluctuation
%   d2gam2, second derivative of structure factor around peak

% results to return
chis=zeros(length(FAV),1);     % spinodal
ks=zeros(length(FAV),1);       % critical wavelength of density fluctuations
d2gam2=zeros(length(FAV),1);   % inverse of susceptibility
NRR=50;

for ii=1:length(FAV)
    fprintf('Step 1: Calculating spinodal at N=%.2e,FA=%.2f\n',N,FAV(ii))
    FA=FAV(ii);

    % find kstar
    G=@(k) gamma2(N,FA,k,0);

    % initial guesses of kstar
    if N<=1e-2
      R2 = r2(NRR);
      k0=-1e-2/(sqrt(R2));
      kf=1e1/(sqrt(R2));
    else
      R2 = r2(N);
      k0=-1e-2/(sqrt(R2));
      kf=1e1/(sqrt(R2));	
    end
    ks(ii) = fminbnd(G,k0,kf);
    chis(ii) = 0.5*G(ks(ii));

    % find susceptibility (second der. of gamma2 w/ k at kstar)
    dks = 1/sqrt(R2)*1e-4;

    if ks(ii)>1e-1  % central differences
        d2gam2(ii) = (G(ks(ii)+dks)-2*G(ks(ii))+G(ks(ii)-dks))/(dks^2);
    else  % forward differences
        d2gam2(ii) = (G(ks(ii)+2*dks)-2*G(ks(ii)+dks)+G(ks(ii)))/(dks^2);
    end
end
