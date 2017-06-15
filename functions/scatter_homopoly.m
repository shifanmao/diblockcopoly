function S = scatter_homopoly(N, KV, CHI, PHIP)
    %%%%%%%%% default inputs for homopolymers in solvent %%%%%%%%%
    FA = 1.0;
    CHIMAT = [CHI, 0;0, 0];
    
    %%%%%%%%% calculate structure factor %%%%%%%%%
    [EIG,~]=gamma2_solvent(N,FA,KV,CHIMAT,PHIP);
    S = 1./EIG(:,1);
end