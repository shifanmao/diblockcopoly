%% calculate mean-field spinodal and critical wavelength at FA=0.5
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

clear;
% Figure 3: mean-field spinodal and critical wavelength at FA=0.5
% NV=logspace(0,3,16);  % number of statistical steps of total chain

V = 0.1;
DSIM = 5;
GV = [32,64];
EPSV = [0.05,0.1,0.5,1.0];
results = [];

for G = GV
    for EPS = EPSV;
        N = EPS*G;

        [chis,ks,d2gam2]=spinodal(N,0.5);

        CHIS = chis*N;
        D = 2*pi/(ks);

        % Simulation input parameters
        L0 = EPS*DSIM/D;
        LP = L0/(2*EPS);

        % Excluded volume parameter
        R2 = r2(N)*(L0/EPS)^2;
        C2 = R2^3/power(G*V,2);

        results = [results;G,EPS,2*LP,C2,CHIS];
        fprintf('G = %d, EPS = %.2f, L0 = %.4f,  2LP = %.4f, R2=%.2f, C2 = %.2f, CHIS = %.2f\n\n', G, EPS, L0, 2*LP, R2, C2, CHIS)
    end
end

% save to file
filename = 'sdata/simparams';
dlmwrite(filename,results,'precision','%.6f');


% data = [log10(NV)',(chis.*NV)',ks'];
% dlmwrite('chivals',data,'precision','%.3f')