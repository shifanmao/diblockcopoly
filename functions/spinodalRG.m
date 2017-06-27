function [chit,phase]=spinodalRG(N,alpha,FAV)
% SPINODALRG  Find the renormalized spinodal by FH theory
% Usage: [chit,phase]=spinodalRG(N,alphaV,FA)
% Inputs:
%   N, number of statistical steps of total chain
%   FA, A monomer chemical composition
% Outputs:
%   chit, renormalized spinodal
%   phase, name of first stable phase above renormalized spinodal,
%     e.g. phase=3, HEXAGONAL phase is the first stable phase

% results to return
chit=zeros(length(FAV),1);     % renormalized spinodal
phase=zeros(length(FAV),1);     % renormalized spinodal

% find spinodal
[chisV,ksV,d2gam2V]=spinodal(N,FAV);

% calculate vertex functions
NQ=1;  % assume to Q dependence
[gam3,gam4]=calcgamma(N,FAV,NQ);
gam3V=real(gam3);
gam4V=real(gam4(:,1));

for ii = 1:length(FAV)
    fprintf('Step 3: Calculating renormalized spinodal at N=%.2e,alpha=%.2f,FA=%.2f\n',N,alpha,FAV(ii))
    
    FA = FAV(ii);
    d2gam2 = d2gam2V(ii);chis = chisV(ii);ks = ksV(ii);
    gam3 = gam3V(ii);gam4 = gam4V(ii);
    if abs(FA-0.5)<1e-2
        chit(ii)=spinodalfh(N,alpha,d2gam2,gam3,gam4,ks,chis,1);
        phase(ii)=1;
    else
        % find renormalized spinodal of each phase
        chi1=spinodalfh(N,alpha,d2gam2,gam3,gam4,ks,chis,1);
        chi3=spinodalfh(N,alpha,d2gam2,gam3,gam4,ks,chis,3);
        chi6=spinodalfh(N,alpha,d2gam2,gam3,gam4,ks,chis,6);

        % find renormalized spinodal
        chiall=[chi1,chi3,chi6];
        phases=[1,3,6];
        [chi,order]=sort(chiall);
        chit(ii)=chi(1);
        phase(ii)=phases(order(1));
    end
end
end

function chit=spinodalfh(N,alpha,d2gam2,gam3,gam4,ks,chis,n)
% calculate renormalized spinodal from FH theory    
    % calculate constant (estimate local second-order derivative)
    c = 0.5*d2gam2;
    Gi = power(alpha,3)*4*pi*sqrt(c)/(ks^2*sqrt(N));
    miu=N*gam3/power(c,3/2);
    lam=N*gam4(1)/power(c,2);
    
    %%%%% LAM/HEX/BCC phase %%%%%
    if n==1
        theta=0*miu;
        eta=-lam/2;
    elseif n==3
        theta=-miu;
        eta=-lam/2;
    elseif n==6
        theta=-2*miu;
        eta=3*lam/2;
    end

    % start solving self-consistent equations
    
    % solver option
    tollev = 1e-16;
    options = optimset('Display','off',...
        'TolX',tollev,'TolFun',tollev,'MaxFunEvals',1/tollev,'MaxIter',1/tollev);

    % initial guesses, lower bounds and upper bounds of solutions
    x0 = [-1e1,1,1,1]/N;
    bnd = 1e5;
    lb=[-bnd,0,0,0];
    ub=[0,bnd,bnd,bnd];
    
    % find corrected spinodal
    % use Hartree approximation from Fredrickson-Helfand theory
    % to evalute free energy of diblock copolymers

    % unknowns x(1)=tau, x(2)=r0, x(3)=r, x(4)=a
    % Equation 1: r0 = tau+sqrt(c)/Gi*power(r0,-1/2)
    % Equation 2: r  = tau+sqrt(c)/Gi*power(r,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               sqrt(c)/Gi*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-sqrt(c)/Gi*lam*power(x(2),-1/2),...
              x(3)-x(1)-sqrt(c)/Gi*lam*power(x(3),-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+sqrt(c)/Gi*(sqrt(x(3))-sqrt(x(2)))-2*n/3*theta*x(4)^3+n/2*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chit=chis-c*tau/(2*N);
end