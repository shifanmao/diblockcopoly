function [chis,chit,phase,chi13,chi36]=plotphaseRG(N,C,FAV)
% PLOTPHASERG:: Calculates diblock copolymer phase diagram
% with Fredrickson-Helfand density fluctuation correction
% Usage :: [chis,chit,phase,chi13,chi36]=plotphaseRG(N,C,FAV)
% Inputs ::
%    FAV, fraction of A-type monomers

% results to return
NFA=length(FAV);
chi13=zeros(NFA,1);      % spinodal (corrected)
chi36=zeros(NFA,1);      % spinodal (corrected)

% calculate mean-field solution
[chis,ks,d2gam2]=spinodal(N,FAV);

% find renormalized spinodal
[chit,phase]=spinodalRG(N,C,FAV);

% calculate vertex functions
NQ=1;  % assume to Q dependence
[gam3,gam4]=calcgamma(N,FAV,NQ);
gam3=real(gam3);
gam4=real(gam4(:,1));

for ii=1:NFA
    fprintf('Step 4: Calculating renormalized OOT phase diag. at N=%.2e,C=%.2e,FA=%.2f\n',N,C,FAV(ii))
    
    % calculate constant (estimate local second-order derivative)
    alpha=power(d2gam2(ii)/2*N/r2(N),1/2);
    d=r2(N)*ks(ii)^2/(4*pi);
    miu=N*gam3(ii)/power(alpha,3);
    lam=N*gam4(ii)/power(alpha,4);

    %%%%% LAM/HEX phase %%%%%
    if (phase(ii)==6 || phase(ii)==3)
        chi13(ii)=chioot(chis(ii),alpha,d,N,C,miu,lam,1,3);
    else
        chi13(ii)=0;
    end
    
    %%%%% HEX/BCC phase %%%%%
    if (phase(ii)==6)    
        chi36(ii)=chioot(chis(ii),alpha,d,N,C,miu,lam,3,6);
    else
        chi36(ii)=0;
    end
end

% make a phase diagram
ind13=find((phase==3 | phase==6) & chi13*N-chit*N>1e-2);
ind36=find((phase==6) & (chi36*N-chit*N>1e-2) & (chi36*N-chi13*N>1e-2));

figure;hold;set(gca,'fontsize',18)
plot(FAV,chit*N,'k-','linewidth',2)
plot(1-FAV,chit*N,'k-','linewidth',2)
FA05 = find(abs(FAV-0.5)<1e-2);
plot(FAV,chis*N,'k--','linewidth',2)
plot(1-FAV,chis*N,'k--','linewidth',2)
if ~isempty(FA05)
    col = 'k';    
    plot(0.5,chis(FA05)*N,'o','color',col,...
    'MarkerSize',8,'MarkerFaceColor',col)

    plot(0.5,chit(FA05)*N,'s','color',col,...
    'MarkerSize',8,'MarkerFaceColor',col)
end

plot(FAV(ind13),chi13(ind13)*N,'r','linewidth',2)
plot(FAV(ind36),chi36(ind36)*N,'b','linewidth',2)
plot(1-FAV(ind13),chi13(ind13)*N,'r','linewidth',2)
plot(1-FAV(ind36),chi36(ind36)*N,'b','linewidth',2)

xlabel('f_A');ylabel('\chi N');box on
xlim([FAV(1),1-FAV(1)]);ylim([5,20])
end

function chi13=chioot(chis,alpha,d,N,C,miu,lam,n1,n2)
    % solver options
    options = optimset('Display','off',...
        'TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',1e14,'MaxIter',1e14);

    % initial guesses, lower bounds and upper bounds of solutions
    x02=[-1e2,1,1,1,1,1,1];
    lb2=[-1e5,0,0,0,0,0,0];
    ub2=[0,1e5,1e5,1e5,1e5,1e5,1e5];
    
    if ((n1==1 && n2==3) || (n1==3 && n2==1))
        % compare LAM/HEX phases
        theta1=0*miu;
        eta1=-lam/2;
        theta2=-miu;
        eta2=-lam/2;
    elseif ((n1==3 && n2==6) || (n1==6 && n2==3))
        % compare HEX/BCC phases
        theta1=-miu;
        eta1=-lam/2;
        theta2=-2*miu;
        eta2=3*lam/2;
    else
        error('Order-order transition between phases not defined')
    end

    % find corrected spinodal
    % use Hartree approximation from Fredrickson-Helfand theory
    % to evalute free energy of diblock copolymers

    % unknowns x(1)=tau, x(2)=r0_1(LAM), x(3)=r_1(LAM), x(4)=a_1(LAM)
    % unknowns           x(5)=r0_3(HEX), x(6)=r_3(HEX), x(7)=a_3(HEX)
    % Equation 1: r0_1 = tau+d*lam*power(r0_1*Nbar,-1/2) (LAM)
    % Equation 2: r0_3 = tau+d*lam*power(r0_3*Nbar,-1/2) (HEX)
    % Equation 3: r_1 = tau+d*lam*power(r_1*Nbar,-1/2)+1*lam*a_1^2 (LAM)
    % Equation 4: r_3 = tau+d*lam*power(r_3*Nbar,-1/2)+3*lam*a_3^2 (HEX)
    % Equation 5: eta_1*a_1^2-theta_1*a+r_1 = 0
    % Equation 6: eta_3*a_3^2-theta_3*a+r_3 = 0
    % Equation 7: phi_1 = phi_3 = (1/2/lam)*(r_1^2-r0_1^2)+...
    %               d*power(Nbar,-1/2)*(r_1^0.5-r0_1^0.5)...
    %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4 - ...
    %                             (1/2/lam)*(r_1^3-r0_1^3)+...
    %               d*power(Nbar,-1/2)*(r_1^0.5-r0_1^0.5)...
    %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4

    F = @(x) [x(2)-x(1)-d/C*lam*power(x(2),-1/2),...
              x(5)-x(1)-d/C*lam*power(x(5),-1/2),...
              x(3)-x(1)-d/C*lam*power(x(3),-1/2)-n1*lam*x(4)^2,...
              x(6)-x(1)-d/C*lam*power(x(6),-1/2)-n2*lam*x(7)^2,...
              eta1*x(4)^2-theta1*x(4)+x(3),...
              eta2*x(7)^2-theta2*x(7)+x(6),...
              ((1/2/lam)*(x(3)^2-x(2)^2)+d/C*(sqrt(x(3))-sqrt(x(2)))-...
                2*n1/3*theta1*x(4)^3+(1/2)*n1*eta1*x(4)^4)-...
              ((1/2/lam)*(x(6)^2-x(5)^2)+d/C*(sqrt(x(6))-sqrt(x(5)))-...
                2*n2/3*theta2*x(7)^3+(1/2)*n2*eta2*x(7)^4)];
    [x,~] = lsqnonlin(F,x02,lb2,ub2,options);
    tau=x(1);
    chi13=chis-alpha^2*tau/(2*N);
end