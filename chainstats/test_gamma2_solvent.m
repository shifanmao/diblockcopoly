clear;

gam2inv = []
eigv = []
CHI = [0,0;0,0];
FAV = linspace(0,1,50);

for FA = FAV
    [EIGS, EIGV] = gamma2_solvent(1,FA,10,CHI,.5);
    gam2inv = [gam2inv, 1/EIGS(1)];
    eigv = [eigv, EIGV(1)];
end

figure;plot(FAV, gam2inv, 'ko-')