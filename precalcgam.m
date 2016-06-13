% Run and save data
clear;
addpath('functions')
addpath('chainstats')
addpath('misc')

FAV = linspace(0.1,0.5,41);  % invariant degree of polymerization
NQ=51;  % number of wavevector sets in calculating GAM4
NV=logspace(-1,4,21);

% write to file
filename='data/newgamdata';
%format  ='%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n';
format  = strcat(repmat('%.2f, ',1,2+NQ),' %.2f\n')
if ~exist(filename,'file')
    outfile = fopen(filename, 'wt');
    for N=NV
        for FA=FAV
            [GAM3,GAM4]=calcgamma(N,FA,NQ);
            fprintf(outfile,format,N,FA,GAM3*N,GAM4(1:NQ)*N);
        end
    end
    fclose(outfile);
end
