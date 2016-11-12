%% calculate mean-field spinodal and critical wavelength at FA=0.5
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

clear;

%%%%%%%%%%%%%%%%%%%%%%%% METHOD 1 :: KEEP R END-TO-END or C CONST. %%%%%%%%%%%%%%%%%%%%%%%%

V = 0.1;
C2V = [1e3,1e4];
GV = [32];
NV = [.05,.1,.2,.5,1,2,5]*32;
results = [];

for C2target = C2V
for G = GV
    for N = NV
        EPS = N/G;
        
        Ree = (sqrt(C2target)*G*0.1)^(1/3);
        L0 = Ree*EPS*((-0.5+0.5*exp(-2*EPS*G)+EPS*G)^(-0.5));

        [chis,ks,d2gam2]=spinodal(N,0.5);

        CHIS = chis*N;

        % Simulation input parameters
        LP = L0/(2*EPS);
        D = (2*pi/ks)*2*LP;
        
        % Excluded volume parameter
        R2 = r2(N)*(L0/EPS)^2;
        C2 = R2^3/power(G*V,2);
        LoverV = power(G*L0,3)/(G*V);
        C22 = r2(N)^3/N^6*LoverV^2; % alternative way to calculate C2 or C squared
        
        results = [results;N,C2,G,2*LP,L0,CHIS];
        fprintf('G=%d, EPS=%.2f, D=%.4f, L0=%.4f, 2LP=%.4f, R2=%.2f, C2=%.2f, LoverV=%.2f, CHIS = %.2f\n\n', ...
                 G,    EPS,      D,      L0,      2*LP,     R2,      C2,      LoverV,      CHIS)
    end
    fprintf('\n')
end
end

% save to file
filename = 'simparams';
dlmwrite(filename,results,'precision','%.6f','delimiter','\t');


% 
% %%%%%%%%%%%%%%%%%%%%%%%% METHOD 2 :: KEEP L0 CONST. %%%%%%%%%%%%%%%%%%%%%%%%
% V = 0.1;
% % DSIM = 5;
% L0tot = 15;
% GV = [32,64];
% EPSV = [0.05,0.1,0.5,1.0];
% results = [];
% 
% for G = GV
%     for EPS = EPSV;
%         N = EPS*G;
% 
%         [chis,ks,d2gam2]=spinodal(N,0.5);
% 
%         CHIS = chis*N;
% 
%         % Simulation input parameters
%         L0 = L0tot/G;
%         LP = L0/(2*EPS);
%         D = (2*pi/ks)*2*LP;
%         
%         % Excluded volume parameter
%         R2 = r2(N)*(L0/EPS)^2;
%         C2 = R2^3/power(G*V,2);
%         LoverV = power(G*L0,3)/(G*V);
%         C22 = r2(N)^3/N^6*LoverV^2; % alternative way to calculate C2 or C squared
%         
%         results = [results;G,EPS,2*LP,C2,CHIS];
%         fprintf('G=%d, EPS=%.2f, D=%.4f, L0=%.4f, 2LP=%.4f, R2=%.2f, C2=%.2f, LoverV=%.2f, CHIS = %.2f\n\n', ...
%                  G,    EPS,      D,      L0,      2*LP,     R2,      C2,      LoverV,      CHIS)
%     end
% end
% 
% % save to file
% filename = 'simparams';
% dlmwrite(filename,results,'precision','%.6f');
% 
% 
% % data = [log10(NV)',(chis.*NV)',ks'];
% % dlmwrite('chivals',data,'precision','%.3f')