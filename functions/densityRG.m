function [k,Smf,Sfh,chit,chis,CHIV,sinvmf,sinvfh]=densityRG(N,alpha,FA,CHIV,PLOTDENSITY,PLOTRG,NCHI)
% function [chis,chit,CHIV,Smf,Sfh,sinvmf,sinvfh]=densityRG(N,C,FA,CHIV,PLOTDENSITY,PLOTRG,NCHI)
% Plot density-density correlations during phase transition
% according to Mean-field theory and FH theory
% Usage :: [chis,chit]=densityRG(N,Nbar,FA)
% Inputs ::
%   N, number of Kuhn steps of whole chain
%   Nbar, invariant degree of polymerization
%   FA, fraction of A monomers
%   CHIV, range of CHI values

% find mean-field spinodal
[chis,ks,d2gamma2]=spinodal(N,FA);

if nargin == 3
    CHIV=(0:0.2:0.8)*chis;
    PLOTRG = 1;
    PLOTDENSITY = 1;
    NCHI = 500;    
elseif nargin == 4
    PLOTRG = 1;
    PLOTDENSITY = 1;
    NCHI = 500;
end

% find quartic order paramter
% gamma4 at angle phi=pi (assume no q dependence)
NQ=1;
[~,gam4]=calcgamma(N,FA,NQ);
gam4=real(gam4(end,1));

% PLOT1: DENSITY CORRELATION  %%
% wavevectors
kmin=0.01;kmax=5;
k=linspace(kmin,kmax,200);

% Flory-Huggins parameter
if (PLOTDENSITY)
    figure;hold;set(gca,'fontsize',20)
    CHIV=(0:0.2:0.8)*chis;
    for ii = 1:length(CHIV)
        CHI=CHIV(ii);
        if length(CHIV)>1
            col = (ii-1)/(length(CHIV)-1);
        else
            col = 0;
        end
        
        Smf=densitymf(N,FA,k,CHI);
        plot(k,Smf,'color',[col 0 1-col],'linestyle','--','linewidth',2);
        plot(ks,densitymf(N,FA,ks,CHI),'o',...
            'MarkerSize',8,'MarkerEdgecolor',[col 0 1-col]);

        % plot RG results
        Sfh=densityfh(N,alpha,FA,k,ks,CHI,d2gamma2,gam4);
        plot(k,Sfh,'color',[col 0 1-col],'linestyle','-','linewidth',2);
        plot(ks,densityfh(N,alpha,FA,ks,ks,CHI,d2gamma2,gam4),'o',...
            'MarkerSize',8,'MarkerFacecolor',[col 0 1-col],'MarkerEdgecolor',[col 0 1-col]);
    end
    xlabel('2l_Pq');
    ylabel('$\left<\psi(\vec{q}^*)\psi(-\vec{q}^*)\right>/Nv$','Interpreter','latex')
	box on
end

% PLOT2: CRITICAL MODE vs CHI
if (PLOTRG)
    % Flory-Huggins parameter
    CHIV=linspace(0,4,NCHI);

    Smf = zeros(length(CHIV),1);
    Sfh = zeros(length(CHIV),1);
    for ii = 1:length(CHIV)
        CHI = CHIV(ii)*chis;
        Smf(ii)=densitymf(N,FA,ks,CHI);
        Sfh(ii)=densityfh(N,alpha,FA,ks,ks,CHI,d2gamma2,gam4);
    end

    figure;hold;set(gca,'fontsize',20)
    col = 'k';
    plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col);
    plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col);

    xlim([1,17]);ylim([0,20]);box on
    xlabel('\chi N');
    ylabel('$N\left<\psi(\vec{q}^*)\psi(-\vec{q}^*)\right>^{-1}$','Interpreter','latex')

    sinvmf=densitymf(N,FA,ks,chis);
    plot(chis*N,1./sinvmf,'o','color',col,...
        'MarkerSize',15,'MarkerFaceColor',col);
    % find renormalized spinodal
    chit=spinodalRG(N,alpha,FA);
    sinvfh=densityfh(N,alpha,FA,ks,ks,chit,d2gamma2,gam4);
        plot(chit*N,1./sinvfh,'s','color',col,...
        'MarkerSize',15,'MarkerFaceColor',col);
end

end

function Smf=densitymf(N,FA,k,CHI)
% Mean-field theory (Leibler)
Gmf=gamma2(N,FA,k,CHI);
Smf=1./(N*Gmf);  % Density-density correlation
end

function Sfh=densityfh(N,alpha,FA,k,ks,CHI,d2gam2,gam4)
% Calculate a few constants gamma2
gam2 = gamma2(N,FA,ks,CHI);

% solve self-consistent equations
c = 0.5*d2gam2;
Gi = power(alpha,3)*4*pi*sqrt(c)/(ks^2*sqrt(N));
tau = N*gam2/c;
lam = N*gam4/c^2;
root = roots([1,0,-tau,-sqrt(c)*lam/Gi]);
r = power(root(root>0),2);

Gfh = c*(r+N*(k-ks).^2);
Sfh=1./(Gfh);  % Density-density correlation
end