function val=gamma3(N,FA,k)
%% This function calculates the cubic order coefficient
% in free energy expansion of copolymer melts
% Usage: G=gamma3(N,FA,k)
% Return:
%    val, Cubic order expansion coeffiecient of free energy
% Parameters:
%    N, number of Kuhn steps per monomer
%    FA, fraction of A type monomers
%    k, magnitude of wavevector in unit of 1/contour length
%       k input can be a vector
% Example:
%     N=1000;
%     k=sqrt(logspace(-1,3,50)/N);
%     FA=0.1;
%     G=gamma3(N,FA,k);
%     figure;hold;set(gca,'fontsize',15);leg=[];
%     plot(k.^2.*NM,G);
%     xlabel('(kR)^2');ylabel('\Gamma_3')
% Shifan Mao 06/10/15

%result to return
val=zeros(length(k),1);

%combination matrix
M=combinator(2,3,'p','r');

for j=1:length(k)
    
    %wave vectors
    Q1=k(j)*[1,0,0];
    Q2=transpose(rotz(pi*2/3)*Q1(1:3)');
    Q3=-Q1-Q2;    

    if N>=1e4  % Gaussian chain limit
        s3 = s3gc(N,FA,Q1,Q2,Q3);
    elseif N<=1e-4  % Rigid rod limit
        s3 = s3rr(N,FA,Q1,Q2,Q3);
    else
        s3 = s3wlc(N,FA,Q1,Q2,Q3);
    end
    s2inv=s2inverse(N,FA,k(j));

    % use Leibler's formula:
    for I = 1:length(M)
        val(j) = val(j) - real(...
                s3(M(I,1),M(I,2),M(I,3))*...
                (s2inv(M(I,1),1)-s2inv(M(I,1),2))*...
                (s2inv(M(I,2),1)-s2inv(M(I,2),2))*...
                (s2inv(M(I,3),1)-s2inv(M(I,3),2)));
    end
end
val=val*power(N,2);
end