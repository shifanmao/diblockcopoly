% use momentum shell integration to provide recursive
% one-loop corrections of correlation functions
function r=shellrg
addpath('functions/')
addpath('chainstats/')
addpath('chainstats/eigcalc/')
addpath('misc/')
clear;

% RG parameters
b = 1.001; % resize scale
Nb = 4000;  % number resizing steps

% Polymer parameters
N = 1e4; % number of Kuhn steps per total chain
FA = 0.50; % fraction of A-type monomers
C = 1e2; % dimensionless volume parameter
R = sqrt(r2(N)); % end-to-end separation

% calculate density-density correlation at peak
[chis,ks,d2gam2]=spinodal(N,FA);
CHI=0.1*chis;
alpha = power(d2gam2/2*N/r2(N),1/2);
tau = (-2*CHI*N+2*chis*N)/alpha^2;

% Calculate mean-field quadratic correlations
lambda0 = 5*ks; % initial cutoff of k
[k,lambda,g2]=initgam2(N,FA,CHI,lambda0,ks,R,b,Nb);
figure;hold;
plot(k(:,1),g2(:,1))
plot(k(:,2),g2(:,2))

% Calculate mean-field vertices
[gam3,gam4]=calcgamma(N,FA,1);

% Start recursive calculation
r = zeros(Nb,1);r(1) = tau;
for ii = 2:Nb
    d = lambda(ii-1)^2*r2(N)/4/pi;
    
    % calculate numerical integral
    int = (lambda(ii-1)-lambda(ii))*(g2(ii-1,1)+g2(ii-1,2));
    r(ii) = b^2*(r(ii-1) + d/C*R*int*N*gam4);
end

end

function [k,lambda,g2]=initgam2(N,FA,CHI,lambda0,ks,R,b,Nb)
    % results to return
    k = zeros(Nb,2);
    g2 = zeros(Nb,2);
    lambda = zeros(Nb,1);

    % initialize wavevectors and quadratic correlation
    lambda(1) = lambda0;

    % lower k cutoff
    ii=1;
    k(ii,1) = ks-lambda(ii);
    if k(ii,1)*R<1e-2
        g2(ii,1) = 0;
    else
        g2(ii,1) = 1/gamma2(N,FA,k(ii,1),CHI);
    end

    % higher k cutoff
    k(ii,2) = ks+lambda(ii);
    g2(ii,2) = 1/gamma2(N,FA,k(ii,2),CHI);

    for ii = 2:Nb
        lambda(ii) = lambda(ii-1)/b;

        % lower k cutoff
        k(ii,1) = ks-lambda(ii);
        if k(ii,1)*R<1e-2
            g2(ii,1) = 0;
        else
            g2(ii,1) = 1/gamma2(N,FA,k(ii,1),CHI);
        end

        % higher k cutoff
        k(ii,2) = ks+lambda(ii);
        g2(ii,2) = 1/gamma2(N,FA,k(ii,2),CHI);
    end
end