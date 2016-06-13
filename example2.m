clear;

% parameters
FA = 0.5;
NV = logspace(0,4,11);
NbarV = logspace(2,10,21);

% results to return
chi1 = zeros(length(NV),length(NbarV));

iN=1;
for N=NV
    iNbar=1;
    for Nbar=NbarV
        chi1(iN,iNbar)=spinodalRG(N,Nbar,FA);
        fprintf('N=%.2f, Nbar=%.2e--> chiN=%.2f\n\n',N,Nbar,chi1(iN,iNbar)*N)
        iNbar = iNbar+1;
    end
    iN = iN+1;
end

figure;hold
iN=1;
for N=NV
    col=(iN-1)/(length(NV)-1);
    plot(NbarV,chi1(iN,:)*N,'-','color',[col 0 1-col]);
    iN=iN+1;
end
plot([NbarV(1),NbarV(end)],[10.495,10.495],'k--')
plot(logspace(log10(NbarV(1)),log10(NbarV(end)),100),...
    10.495+41.022*power(logspace(log10(NbarV(1)),log10(NbarV(end)),100),-1/3),'k-')
xlim([NbarV(1),NbarV(end)])
set(gca,'xscale','log')

figure;hold
iNbar=1;
for Nbar=NbarV
    col=(iNbar-1)/(length(NbarV)-1);
    plot(NV,chi1(:,iNbar).*NV','-','color',[col 0 1-col]);
    iNbar=iNbar+1;
end
set(gca,'xscale','log')