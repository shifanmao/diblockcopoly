clear;
addpath('../functions/')
data = load('../data/gamdata');
cd ../

NV = unique(data(:,1))';

% NV=logspace(-1,4,21);
% FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
FAV=0.5;

[chis,ks,d2gam2]=spinodal(NV,0.5);
alpha=power(d2gam2/2.*NV./r2(NV),1/2);

gam3 = zeros(length(FAV),1);
gam4 = zeros(length(FAV),1);
NQ=1;

ii=1;
for N = NV
    [gam3(ii),gam4(ii)]=calcgamma(N,FAV,NQ);
    ii = ii+1;
end