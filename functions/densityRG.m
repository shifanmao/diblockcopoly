function [chis,chit,CHIV,Smf,Sfh,sinvmf,sinvfh]=densityRG(N,C,FA)
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

% PLOT1: DENSITY CORRELATION  %%
% wavevectors
kmin=0.1;kmax=20;
k=power(linspace(kmin,kmax,200),1)/sqrt(r2(N));

% Flory-Huggins parameter
CHIV=0:0.2:0.8;

figure;hold;set(gca,'fontsize',20)
p1=[];p2=[];
for ii = 1:length(CHIV)
    CHI=chis*CHIV(ii);
    if length(CHIV)>1
        col = (ii-1)/(length(CHIV)-1);
    else
        col = 0;
    end

    if CHI<=0.9*chis
        % plot mean-field results
        Smf=densitymf(N,FA,k,CHI);
        p1(ii)=plot(k*sqrt(r2(N)),Smf./N,'color',[col 0 1-col],'linestyle','--','linewidth',2);
        plot(ks*sqrt(r2(N)),densitymf(N,FA,ks,CHI)/N,'o',...
            'MarkerSize',8,'MarkerEdgecolor',[col 0 1-col]);
    end
    
    % plot RG results
    Sfh=densityfh(N,C,FA,k,ks,CHI,d2gamma2);
    p2(ii)=plot(k*sqrt(r2(N)),Sfh./N,'color',[col 0 1-col],'linestyle','-','linewidth',2);
    plot(ks*sqrt(r2(N)),densityfh(N,C,FA,ks,ks,CHI,d2gamma2)/N,'o',...
        'MarkerSize',8,'MarkerFacecolor',[col 0 1-col],'MarkerEdgecolor',[col 0 1-col]);
end
xlim([kmin,kmax]);box on
xlabel('qR');ylabel('<\psi^2>/N')
% legend(p2,{'\chi=0','\chi=0.2\chi_S^{MF}',...
%     '\chi=0.4\chi_S^{MF}','\chi=0.6\chi_S^{MF}','\chi=0.8\chi_S^{MF}'},'location','northeast')

% PLOT2: CRITICAL MODE vs CHI
% Flory-Huggins parameter
CHIV=linspace(0,4,200);

Smf = zeros(length(CHIV),1);
Sfh = zeros(length(CHIV),1);
for ii = 1:length(CHIV)
    CHI = CHIV(ii)*chis;
    Smf(ii)=densitymf(N,FA,ks,CHI);
    Sfh(ii)=densityfh(N,C,FA,ks,ks,CHI,d2gamma2);
end

figure;hold;set(gca,'fontsize',20)
col = 'k';
plot(CHIV*chis*N,1./Smf,'--','linewidth',2,'color',col);
plot(CHIV*chis*N,1./Sfh,'-','linewidth',2,'color',col);
xlim([1,17]);ylim([0,20]);box on
% xlabel('\chi N');ylabel('1/<\psi^2(q^*)>')
xlabel('\chi N');ylabel('$N<\tilde{\psi}(q^*)\tilde{\psi}(-q^*)>^{-1}$','Interpreter','latex')

% find renormalized spinodal
chit=spinodalRG(N,C,FA);
sinvmf=densitymf(N,FA,ks,chis);
sinvfh=densityfh(N,C,FA,ks,ks,chit,d2gamma2);
plot(chis*N,1./sinvmf,'o','color',col,...
    'MarkerSize',8,'MarkerFaceColor',col);
plot(chit*N,1./sinvfh,'s','color',col,...
    'MarkerSize',8,'MarkerFaceColor',col);
end

function Smf=densitymf(N,FA,k,CHI)
% Mean-field theory (Leibler)
Gmf=gamma2(N,FA,k,CHI);
Smf=1./(N*Gmf);  % Density-density correlation
end

function Sfh=densityfh(N,C,FA,k,ks,CHI,d2gamma2)
% RG with no q dependence (F-H)
% self-consistent equations:
%  unknowns x(1)=r, x(2)=alpha
%  r = N*gamma2+pref*N*gamma4
%  2*alpha = d^2(gamma2)/dQ^2 + pref*d^2(gamma4)/dQ^2
%       , where pref = N*Qs/(4*pi*power(Nbar*r*alpha,1/2))

% Calculate a few constants gamma2
gam2 = gamma2(N,FA,ks,CHI);

% gamma4 at angle phi=pi (assume no q dependence)
NQ=1;
[~,gam4]=calcgamma(N,FA,NQ);
gam4=real(gam4(end,1));

% solve self-consistent equations
alpha=power(d2gamma2/2*N/r2(N),1/2);
d=r2(N)*ks^2/(4*pi);
root = roots([1,0,-N*gam2/(alpha^2),-(d/C)*N*gam4/(alpha^4)]);
r = power(root(root>0),2);

Gfh = alpha^2*(r + r2(N)*(k-ks).^2);
Sfh=1./(Gfh);  % Density-density correlation
end

% Reserved for future work
% %% RG with no q dependence (F-H)
% CHAIN=1;
% ORDEig=4;NumLayer=500;
% % self-consistent equations:
% %  unknowns x(1)=r, x(2)=alpha
% %  r = N*gamma2+pref*N*gamma4
% %  2*alpha = d^2(gamma2)/dQ^2 + pref*d^2(gamma4)/dQ^2
% %       , where pref = N*Qs/(4*pi*power(Nbar*r*alpha,1/2))
% 
% % Calculate a few constants
% % gamma2
% gam2 = gamma2(CHAIN,N,NM,FA,Qs,0,ORDEig,NumLayer);
% % second-order derivative of gamma2 near Qs
% dQ=Qs*1e-3;
% d2gamma2=(gamma2(CHAIN,N,NM,FA,Qs+dQ,0,ORDEig,NumLayer)-...
%             2*gamma2(CHAIN,N,NM,FA,Qs,0,ORDEig,NumLayer)+...
%             gamma2(CHAIN,N,NM,FA,Qs-dQ,0,ORDEig,NumLayer))/(dQ^2);
% % gamma4 at angle phi=pi (assume no q dependence)
% Q1=[1,0,0]';Q2=rotz(pi)*Q1;Q3=-Q2;Q4=-Q1;
% gam4 = gamma4(Qs,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer);
% d2gamma4 = (gamma4(Qs+dQ,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer)-...
%             2*gamma4(Qs,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer)+...
%             gamma4(Qs-dQ,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer))/(dQ^2);
%         
% % % initial guesses, lower bounds and upper bounds of solutions
% % unknowns: x = [r,alpha]
% x0 = [1,1];
% lb=[0,0];
% ub=[1e5,1e5];
% options = optimset('Display','off',...
%     'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);
% F = @(x) [x(1)-N*gam2-(N*Qs^2)/(4*pi*power(Nbar*x(1)*x(2),1/2))*N*gam4,...
%           2*x(2)-d2gamma2-(N*Qs^2)/(4*pi*power(Nbar*x(1)*x(2),1/2))*d2gamma4];
% [x,~] = lsqnonlin(F,x0,lb,ub,options);
% r = x(1); alpha = x(2);
% Gfh = r/N + alpha*(Q-Qs).^2;
% chisfh = (r/N)/2*N

% % RG with q dependence (Fredrickson-Barret)
% % self-consistent equations:
% %  unknowns x(1)=r, x(2)=Qstar, x(3)=alpha
% %  r = N*gamma2+pref*N*gamma4
% %  0 = d(gamma2)/dQ + pref*d(gamma4)/dQ
% %  2*alpha = d^2(gamma2)/dQ^2 + pref*d^2(gamma4)/dQ^2
% %
% %  where pref = N*Qstar^2/(4*pi*power(Nbar*r*alpha,1/2))
% r = 0; alpha = 0; Qstar = Qs;
%     
% % initial guesses, lower bounds and upper bounds of solutions
% x0 = [1,Qs,1];
% lb=[0,0,0];
% ub=[1e5,1e5,1e5];
% options = optimset('Display','off',...
%     'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);
% F = @(x) [x(1)-N*gamma2-(N*x(2)^2)/(4*pi*power(Nbar*x(1)*x(3),1/2))*N*gamma4,...
%           dgamma2+(N*x(2)^2)/(4*pi*power(Nbar*x(1)*x(3),1/2))*dgamma4,...
%           x(3)-d2gamma2-(N*x(2)^2)/(4*pi*power(Nbar*x(1)*x(3),1/2))*d2gamma4];
% [x,~] = lsqnonlin(F,x0,lb,ub,options);
% G = r/N + alpha*(Q-Qstar).^2;
