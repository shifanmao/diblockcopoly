clear;

N = 1;
FA = 0.5;

PHIPV = [0.1:0.1:1];
CHIV = linspace(0, .2, 101);
KSV = zeros(length(PHIPV), length(CHIV));

for ii = 1:length(PHIPV)
    PHIP = PHIPV(ii)
    [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N, FA, PHIP);

    for jj = 1:length(CHIV)
        CHI = CHIV(jj);
        CHIMAT = [0,CHIABS;CHIABS,0]*CHI;
        [KS, EIGS, EIGVS] = eigs_wsolvent(N, FA, PHIP, CHIMAT, 1);

        if KS*sqrt(r2(N))<.1+1e-4
            KS = 0;
        end
        KSV(ii,jj) = KS;
    end
end

figure;hold
for ii = 1:length(PHIPV)
    COL = (ii-1) / (length(PHIPV)-1);
    plot(CHIV, KSV(ii,:)*sqrt(r2(N)), 'o-',...
        'color', [COL 0 1-COL], 'linewidth', 2)
end
box on
xlabel('\chi_{AB}/\chi_{AB}^*')
ylabel('Rq^*(\chi_{AB})')
set(gca,'fontsize',18)