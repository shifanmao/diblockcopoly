function [chis,ks,d2gam2]=spinodal(NV,FAV)
% spinodal.m :: This function predicts diblock copolymer phase transition spinodal
% and critical wavelength of quadratic density fluctuation in the well-mixed state
% Usage: [chis,ks,d2gam2]=spinodal(N,FAV)
% Inputs:
%   NV, number of statistical steps of total chain
%   FAV, range of A-type monomer fractions
% Outputs:
%   chis, Flory-Huggins parameter at spinodal
%   ks, critical wavelength of quadratic density fluctuation
%   d2gam2, second derivative of structure factor around peak

% results to return
chis=zeros(length(FAV),length(NV));     % spinodal
ks=zeros(length(FAV),length(NV));        % critical wavelength of density fluctuations
d2gam2=zeros(length(FAV),length(NV));   % inverse of susceptibility

for ii=1:length(FAV)
    FA=FAV(ii);
    for jj=1:length(NV)
        N=NV(jj);
        fprintf('Step 1: Calculating spinodal at N=%.2e,FA=%.2f\n',N,FA)

        % find kstar
        G=@(k) gamma2(N,FA,k,0);

        % initial guesses of kstar
        R2 = r2(N);
        k0=-1e-2/(sqrt(R2));
        kf=1e1/(sqrt(R2));
        ks(ii,jj) = fminbnd(G,k0,kf);

        chis(ii,jj) = 0.5*G(ks(ii,jj));

        % find susceptibility (second der. of gamma2 w/ k at kstar)
        dks = 1/sqrt(R2)*1e-5;

        if ks(ii,jj)>1e-1  % central differences
            d2gam2(ii,jj) = (G(ks(ii,jj)+dks)-2*G(ks(ii,jj))+G(ks(ii,jj)-dks))/(dks^2);
        else  % forward differences
            d2gam2(ii,jj) = (G(ks(ii,jj)+2*dks)-2*G(ks(ii,jj)+dks)+G(ks(ii,jj)))/(dks^2);
        end
    end
end