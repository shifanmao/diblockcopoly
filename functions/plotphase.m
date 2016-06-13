function [chis,chi13,chi36,chi12,chi23,chi26]=plotphase(N,FAV)
% Plots phase mean-field diagram of diblock copolymers
% Usage: plotphase(N,FAV)
% Parameters:
%   N, number of Kuhn steps of total chain
%   FAV, fraction of A-type monomers
% Output:
%   chis, ODT Flory-Huggins parameter
%   chi13, OOT between Lamellar (1) and Cylinder (3) phase
%   chi36, OOT between Cylinder (3) and Body-centered cubic (6) phase
%   , and a plot of phase diagram

[chis,~]=spinodal(N,FAV);
[gam3,gam4]=calcgamma(N,FAV,4);

% make sure imaginary parts are small
gam3= real(gam3);
gam4 = real(gam4);
[chi13,chi36,chi12,chi23,chi26]=phasediag(N,FAV,gam3,gam4);

chi13(chi13<1e-8)=0;
chi36(chi36<1e-8)=0;

% make a phase diagram
figure;set(gca,'fontsize',18);hold
col1='k-';
col2='r-';
col3='b-';

plot(FAV,chis*N,col1,'linewidth',2);
plot(1-FAV,chis*N,col1,'linewidth',1.5)
FA05 = find(abs(FAV-0.5)<1e-2);
if ~isempty(FA05)
    col = 'k';
    plot(0.5,chis(FA05)*N,'o','color',col,...
    'MarkerSize',8,'MarkerFaceColor',col)
end

plot(FAV,(chi13+chis)*N,col2,'linewidth',1.5)
plot(1-FAV,(chi13+chis)*N,col2,'linewidth',1.5)

plot(FAV,(chi36+chis)*N,col3,'linewidth',1.5)
plot(1-FAV,(chi36+chis)*N,col3,'linewidth',1.5)
xlabel('f_A');ylabel('\chi N');box on
xlim([FAV(1),1-FAV(1)]);ylim([5,20]);
end

function [chi13,chi36,chi12,chi23,chi26]=phasediag(N,FAV,gam3,gam4)
% calculates order-order transition (OOT) boundaries
% based on free energy expansion of each phase
    NF = length(FAV);
    chi13=zeros(NF,1);
    chi36=zeros(NF,1);
    chi12=zeros(NF,1);
    chi23=zeros(NF,1);
    chi26=zeros(NF,1);

    % free energy
    for I=1:NF
        FA=FAV(I);
        
        gamma3=gam3(I);
        gamma4=gam4(I,:);
        
        fprintf('Step 2: Calculating OOT phase diag. at FA=%.2f, N=%.2e\n',FA,N)
        F13=@(CHI) abs(lamellar(CHI,gamma4)-hexagonal(CHI,gamma3,gamma4));
        F36=@(CHI) abs(hexagonal(CHI,gamma3,gamma4)-bcc(CHI,gamma3,gamma4));
        F12=@(CHI) abs(lamellar(CHI,gamma4)-tetragonal(CHI,gamma4));
        F23=@(CHI) abs(tetragonal(CHI,gamma4)-hexagonal(CHI,gamma3,gamma4));
        F26=@(CHI) abs(tetragonal(CHI,gamma4)-bcc(CHI,gamma3,gamma4));
        
        options = optimset('Display','off',...
            'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);
        chi13(I) = fminbnd(F13,0,200/N,options);
        chi36(I) = fminbnd(F36,0,200/N,options);
        chi12(I) = fminbnd(F12,0,200/N,options);
        chi23(I) = fminbnd(F23,0,200/N,options);
        chi26(I) = fminbnd(F26,0,200/N,options);
    end
end

function F1=lamellar(CHI,gam4)
    %%% LAM phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B1=0;
    C1=power(2*pi,-9)*(1/4)*gam4(1);

    % free energy %
    FE=@(psi) A*power(psi,2)+B1*power(psi,3)+C1*power(psi,4);
    x1=fminbnd(FE,1e-10,1e5);
    F1=FE(x1);
end
function F2=tetragonal(CHI,gam4)
    %%% TETRA phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B2=0;
    C2=power(2*pi,-9)*(1/8)*(gam4(1)+2*gam4(3));

    % free energy %
    FE=@(psi) A*power(psi,2)+B2*power(psi,3)+C2*power(psi,4);
    x2=fminbnd(FE,1e-10,1e5);
    F2=FE(x2);
end
function F3=hexagonal(CHI,gam3,gam4)
    %%% HEX phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B3=power(2*pi,-6)*(2/3/sqrt(3))*gam3;
    C3=power(2*pi,-9)*(1/12)*(gam4(1)+4*gam4(2));

    % free energy %
    FE=@(psi) A*power(psi,2)+B3*power(psi,3)+C3*power(psi,4);
    x3=fminbnd(FE,1e-10,1e5);
    F3=FE(x3);
end
function F6=bcc(CHI,gam3,gam4)
    %%% BCC phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B6=power(2*pi,-6)*(4/3/sqrt(6))*gam3;
    C6=power(2*pi,-9)*(1/24)*(gam4(1)+8*gam4(2)+2*gam4(3)+4*gam4(4));

    % free energy %
    FE=@(psi) A*power(psi,2)+B6*power(psi,3)+C6*power(psi,4);
    x6=fminbnd(FE,1e-10,1e5);
    F6=FE(x6);
end
