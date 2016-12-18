function val=gamma2(N,FA,k,CHI)
%% This function calculates the quadratic order coefficient
% in free energy expansion of copolymer melts
% Usage: val=gamma2(N,FA,k,CHI)
% Return:
%     val, quadratic order expansion coeffiecient of free energy
% Parameters:
%    k, magnitude of wavevector in unit of 1/contour length
%       k input can be a vector
%    N, number of monomers
%    FA, fraction of A type monomers
%    CHI, chemical incompatibility between A and B monomers, non-
%        dimensionalized by monomer volume v (CHI*v)
% Example:
%    N=1e5;FA=0.5;
%    k = logspace(-2,2,50);G=gamma2(N,FA,k,0)
%    figure;hold;set(gca,'fontsize',18);
%    plot(k,1./G);
%    xlabel('(kR)^2');ylabel('1/\Gamma_2')
% Shifan Mao 06/10/15

%result to return
val=zeros(length(k),1);
NRR=50;

%combination matrix
M=combinator(2,2,'p','r');

%sign indicator
D=[1,-1];

%% evaluate s2inv
for j = 1:length(k)
    s2inv=s2inverse(N,FA,k(j));
    G=0;
    for I = 1:4
        G = G + real(...
                s2inv(M(I,1),M(I,2))*...
                D(M(I,1))*D(M(I,2)));
    end

    if N<=1e-2 
      val(j) = -2*CHI+NRR*G;
    else
      val(j) = -2*CHI+N*G;
    end
end

end
